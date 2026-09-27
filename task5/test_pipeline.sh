#!/bin/bash
# Quick test of main.nf on chromosome 22 + a few thousand reads (catches bugs before the full run)
set -e
source ~/miniforge3/etc/profile.d/conda.sh && conda activate genomics
mkdir -p ~/task5_test/fastq && cd ~/task5_test
R=~/task3/reference
[ -f chr22.fa ]  || samtools faidx $R/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 22 > chr22.fa
[ -f chr22.gtf ] || awk '$1=="22"' $R/Homo_sapiens.GRCh38.91.gtf > chr22.gtf
cp ~/task3/fastq/SGNex_Hct116_directRNA_replicate6_run1.fastq.gz fastq/
zcat ~/task3/fastq/SGNex_A549_directRNA_replicate1_run1.fastq.gz | head -n 80000 | gzip > fastq/SGNex_A549_directRNA_replicate1_run1.fastq.gz
zcat ~/task3/fastq/SGNex_K562_cDNA_replicate1_run3.fastq.gz   | head -n 80000 | gzip > fastq/SGNex_K562_cDNA_replicate1_run3.fastq.gz
ARGS="--reads '$PWD/fastq/*.fastq.gz' --refFa $PWD/chr22.fa --refGtf $PWD/chr22.gtf --threads 2 --qc_min_mapped 0"
eval nextflow run ~/task5/main.nf -ansi-log false --outdir $PWD/results1 $ARGS
eval nextflow run ~/task5/main.nf -ansi-log false --outdir $PWD/results2 $ARGS --annotations false -resume
echo "=== QC summary"; cat results1/qc/qc_summary.tsv
echo "=== outputs"; ls results1/bambu results2/bambu results1/bam
echo TEST_PASSED
