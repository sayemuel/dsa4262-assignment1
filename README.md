# DSA4262 Assignment 1 – long read RNA-Seq

Code for Assignment 1 (Introduction to Genomics). The data is from the SG-NEx project: https://github.com/GoekeLab/sg-nex-data

The Nextflow pipeline in `task5/` is based on the workshop pipeline (`workflow_longReadRNASeq.nf` from GoekeLab/sg-nex-data). I changed it so that it:

- runs all the fastq files at once instead of one sample
- picks the minimap2 settings from the file name (direct RNA: `-ax splice -uf -k14`, cDNA: `-ax splice`)
- sorts and indexes the bam files
- has a QC step (seqkit stats + samtools flagstat) that marks each sample PASS or FAIL, and only PASS samples go to Bambu
- can run Bambu with or without the annotation (`--annotations true/false`)

## Files

- `task3/align.sh` – alignment of the 4 samples with minimap2 + samtools (run on RONIN)
- `task4/run_bambu.R` – Bambu on the 4 bam files, plus the code for the answers to 4.2
- `task5/main.nf` – the Nextflow pipeline
- `task5/test_pipeline.sh` – small test on chr22 with a few thousand reads, to check the pipeline before the full run
- `task5/report_scenario1.html`, `report_scenario2.html` (and timelines) – Nextflow reports for the two runs
- `task5/qc_summary.tsv` – QC results

## Running it

Tools I used (installed in one conda env): nextflow, openjdk 17, minimap2, samtools, seqkit, R with bambu.

```
mamba create -n genomics -c conda-forge -c bioconda nextflow openjdk=17 minimap2 samtools seqkit bioconductor-bambu r-biocmanager
```

Scenario 1 (with annotation):

```
nextflow run task5/main.nf -with-report report_scenario1.html \
  --reads '/path/to/fastq/*.fastq.gz' \
  --refFa /path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --refGtf /path/to/Homo_sapiens.GRCh38.91.gtf \
  --outdir results_scenario1
```

Scenario 2 (no annotation, reusing the cached alignment):

```
nextflow run task5/main.nf -resume -with-report report_scenario2.html \
  --reads '/path/to/fastq/*.fastq.gz' \
  --refFa /path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --refGtf /path/to/Homo_sapiens.GRCh38.91.gtf \
  --outdir results_scenario2 --annotations false
```

Notes:
- put quotes around the `--reads` pattern, otherwise the shell expands it and only one file gets passed in
- the minimap2 index for the whole human genome needs around 20 GB RAM, so the alignment step runs one sample at a time. I ran everything on a RONIN machine with 8 CPUs and 32 GB RAM
- QC thresholds can be changed with `--qc_min_reads`, `--qc_min_meanlen`, `--qc_min_meanq`, `--qc_min_mapped`

## Data

- reads: `s3://sg-nex-data/data/sequencing_data_ont/fastq/` (A549 directRNA rep1 run1, Hct116 directRNA rep6 run1, Hct116 cDNA rep1 run6, K562 cDNA rep1 run3)
- genome: `s3://sg-nex-data/data/annotations/genome_fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`
- annotation: `s3://sg-nex-data/data/annotations/gtf_file/Homo_sapiens.GRCh38.91.gtf`
