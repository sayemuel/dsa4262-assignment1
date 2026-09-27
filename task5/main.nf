#!/usr/bin/env nextflow
// Long-read RNA-Seq pipeline: QC -> minimap2 -> samtools -> Bambu
// Extended from the SG-NEx workshop pipeline (GoekeLab/sg-nex-data, docs/colab/workflow_longReadRNASeq.nf):
//  - runs any number of samples in parallel, with minimap2 flags chosen per protocol (directRNA / cDNA)
//  - sorts and indexes the bam files
//  - adds read-level and alignment-level QC, and only passes samples that pass QC on to Bambu
//  - one switch (--annotations) to run Bambu with or without the reference annotations

params.reads       = "$projectDir/fastq/*.fastq.gz"
params.refFa       = '/path/to/ref.fa'
params.refGtf      = '/path/to/ref.gtf'
params.outdir      = 'results'
params.annotations = true      // scenario 1: true, scenario 2: false
params.threads     = 8
params.bambu_ncore = 1       // Bambu's multi-core mode crashed on the full dataset, so it runs samples one at a time

// QC thresholds (a sample must meet all of them to be used for Bambu)
params.qc_min_reads   = 10000  // enough reads to say anything about expression
params.qc_min_meanlen = 300    // most transcripts are > 300 nt; much shorter reads are mostly fragments
params.qc_min_meanq   = 7      // Nanopore's usual "pass" read quality cut-off
params.qc_min_mapped  = 50     // % of reads that align to the genome
params.qc_filter      = true   // drop failing samples before Bambu


// ---- QC 1: read-level statistics (count, length, N50, quality) ----
process READ_QC {
  tag "$sample"
  publishDir "${params.outdir}/qc", mode: 'copy'
  cpus 2

  input:
    tuple val(sample), path(reads)
  output:
    tuple val(sample), path("${sample}.readstats.tsv")

  script:
  """
  seqkit stats -a -T -j ${task.cpus} $reads > ${sample}.readstats.tsv
  """
}

// ---- alignment: spliced alignment with protocol-specific flags ----
process MINIMAP2_ALIGN {
  tag "$sample"
  cpus params.threads
  maxForks 1              // the full human splice index needs ~22 GB RAM, so align one sample at a time

  input:
    path refFa
    tuple val(sample), path(reads)
  output:
    tuple val(sample), path("${sample}.sam")

  script:
  def flags = sample.contains('directRNA') ? '-ax splice -uf -k14' : '-ax splice'
  """
  minimap2 $flags -t ${task.cpus} $refFa $reads > ${sample}.sam
  """
}

// ---- sam -> sorted, indexed bam ----
process SAM_TO_BAM {
  tag "$sample"
  publishDir "${params.outdir}/bam", mode: 'copy'
  cpus 4

  input:
    tuple val(sample), path(sam)
  output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

  script:
  """
  samtools view -@ ${task.cpus} -b $sam > ${sample}.unsorted.bam
  samtools sort -@ ${task.cpus} ${sample}.unsorted.bam -o ${sample}.bam
  samtools index ${sample}.bam
  rm ${sample}.unsorted.bam
  """
}

// ---- QC 2: alignment statistics ----
process ALIGN_QC {
  tag "$sample"
  publishDir "${params.outdir}/qc", mode: 'copy'
  cpus 2

  input:
    tuple val(sample), path(bam), path(bai)
  output:
    tuple val(sample), path("${sample}.flagstat.txt")

  script:
  """
  samtools flagstat -@ ${task.cpus} $bam > ${sample}.flagstat.txt
  """
}

// ---- QC 3: combine metrics and decide PASS / FAIL ----
process QC_CHECK {
  tag "$sample"

  input:
    tuple val(sample), path(readstats), path(flagstat)
  output:
    tuple val(sample), path("${sample}.qc.tsv"), env('QC_STATUS')

  script:
  """
  # read metrics from seqkit (look columns up by name)
  eval \$(awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++)c[\$i]=i} NR==2{print "NREADS="\$c["num_seqs"]" MEANLEN="\$c["avg_len"]" N50="\$c["N50"]" MEANQ="\$c["AvgQual"]}' $readstats)
  MAPPED=\$(grep "primary mapped" $flagstat | sed -E 's/.*\\(([0-9.]+)%.*/\\1/')

  QC_STATUS=PASS; REASON=""
  awk "BEGIN{exit !(\$NREADS  < ${params.qc_min_reads})}"   && { QC_STATUS=FAIL; REASON="\$REASON low_read_count"; }
  awk "BEGIN{exit !(\$MEANLEN < ${params.qc_min_meanlen})}" && { QC_STATUS=FAIL; REASON="\$REASON short_reads"; }
  awk "BEGIN{exit !(\$MEANQ   < ${params.qc_min_meanq})}"   && { QC_STATUS=FAIL; REASON="\$REASON low_quality"; }
  awk "BEGIN{exit !(\$MAPPED  < ${params.qc_min_mapped})}"  && { QC_STATUS=FAIL; REASON="\$REASON low_mapping_rate"; }

  printf "sample\\treads\\tmean_length\\tN50\\tmean_quality\\tpct_mapped\\tqc_status\\treason\\n" > ${sample}.qc.tsv
  printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" $sample \$NREADS \$MEANLEN \$N50 \$MEANQ \$MAPPED \$QC_STATUS "\${REASON:- -}" >> ${sample}.qc.tsv
  """
}

// ---- QC summary table for all samples ----
process QC_REPORT {
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
    path qc_files
  output:
    path "qc_summary.tsv"

  script:
  """
  head -n 1 \$(ls *.qc.tsv | head -n 1) > qc_summary.tsv
  for f in *.qc.tsv; do tail -n +2 \$f; done | sort >> qc_summary.tsv
  """
}

// ---- transcript discovery and quantification ----
process BAMBU {
  publishDir "${params.outdir}/bambu", mode: 'copy'
  cpus params.threads

  input:
    path refFa
    path refGtf
    path bams
  output:
    path "*.txt"
    path "*.gtf"

  script:
  // command-line values arrive as text ("false" would count as true), so convert explicitly
  def useAnnotations = params.annotations.toString().toBoolean()
  def run = useAnnotations ?
    """se <- bambu(reads = bams, annotations = prepareAnnotations("$refGtf"), genome = "$refFa", ncore = ${params.bambu_ncore})""" :
    """se <- bambu(reads = bams, annotations = NULL, genome = "$refFa", NDR = 1,
                   opt.discovery = list(min.readFractionByEqClass = 0.2), ncore = ${params.bambu_ncore})"""
  """
  #!/usr/bin/env Rscript
  library(bambu)
  bams <- list.files(".", pattern = "\\\\.bam\$")
  print(bams)
  $run
  writeBambuOutput(se, path = "./")
  writeLines(capture.output(show(se)), "bambu_summary.txt")
  """
}

workflow {
  reads = Channel.fromPath(params.reads, checkIfExists: true)
                 .map { f -> tuple(f.name.replaceAll(/\.fastq\.gz$|\.fq\.gz$|\.fastq$/, ''), f) }

  READ_QC(reads)
  MINIMAP2_ALIGN(params.refFa, reads)
  SAM_TO_BAM(MINIMAP2_ALIGN.out)
  ALIGN_QC(SAM_TO_BAM.out)

  QC_CHECK(READ_QC.out.join(ALIGN_QC.out))
  QC_REPORT(QC_CHECK.out.map { s, tsv, status -> tsv }.collect())

  QC_CHECK.out
    .map { s, tsv, status -> tuple(s, status) }
    .subscribe { s, status -> if (status == 'FAIL') log.warn "QC FAIL: $s (see ${params.outdir}/qc/qc_summary.tsv)" }

  bams_for_bambu = SAM_TO_BAM.out
    .join(QC_CHECK.out.map { s, tsv, status -> tuple(s, status) })
    .filter { s, bam, bai, status -> status == 'PASS' || !params.qc_filter }
    .map { s, bam, bai, status -> bam }
    .collect(sort: true)   // fixed order, so -resume does not rerun Bambu just because samples finished in a different order

  BAMBU(params.refFa, params.refGtf, bams_for_bambu)
}
