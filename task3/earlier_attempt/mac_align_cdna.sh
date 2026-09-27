#!/bin/bash
# Task 3 (cDNA samples) on the Mac: minimap2 -> samtools, streamed so no SAM file hits the (nearly full) disk.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh && conda activate genomics

cd ~/genomics_assignment/task3
mkdir -p fastq bam
REF=../reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
T=6   # leave 2 cores free so the Mac stays usable

for s in SGNex_Hct116_cDNA_replicate1_run6 SGNex_K562_cDNA_replicate1_run3; do
  if [ -f bam/$s.bam.bai ]; then echo "=== $s already done, skipping"; continue; fi
  [ -f fastq/$s.fastq.gz ] || aws s3 cp --no-sign-request --only-show-errors \
      s3://sg-nex-data/data/sequencing_data_ont/fastq/$s/$s.fastq.gz fastq/

  echo "=== $(date '+%H:%M:%S')  $s"
  # cDNA protocol: spliced alignment, reads can come from either strand -> no -uf, default k-mer (15)
  # -I 1g --split-prefix: build the genome index in ~1 Gb chunks so it fits in 16 GB RAM
  echo "\$ minimap2 -ax splice -t $T -I 1g --split-prefix bam/tmp_$s $REF fastq/$s.fastq.gz | samtools view -b -o bam/$s.unsorted.bam"
  minimap2 -ax splice -t $T -I 1g --split-prefix bam/tmp_$s $REF fastq/$s.fastq.gz | samtools view -@ 2 -b -o bam/$s.unsorted.bam
  echo "\$ samtools sort bam/$s.unsorted.bam -o bam/$s.bam"
  samtools sort -@ $T bam/$s.unsorted.bam -o bam/$s.bam
  echo "\$ samtools index bam/$s.bam"
  samtools index bam/$s.bam
  rm -f bam/$s.unsorted.bam
  echo "=== $(date '+%H:%M:%S')  $s done"
done
echo "ALL DONE"
ls -lh bam/
