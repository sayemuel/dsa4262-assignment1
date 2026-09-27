#!/bin/bash
# Task 3 on RONIN: align SG-NEx reads with minimap2, convert/sort/index with samtools.
# Usage:  bash ronin_align.sh [sample ...]    (default: all four samples)
# Set AUTO_STOP=1 to shut the machine down when finished (stops billing; disk and files are kept).
set -euo pipefail

WD=~/task3
mkdir -p $WD/reference $WD/fastq $WD/sam $WD/bam $WD/software
cd $WD

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
  SAMPLES=(SGNex_Hct116_directRNA_replicate6_run1 SGNex_A549_directRNA_replicate1_run1
           SGNex_Hct116_cDNA_replicate1_run6 SGNex_K562_cDNA_replicate1_run3)
fi

echo "=== machine: $(nproc) CPUs, $(free -g | awk '/Mem:/{print $2}') GB RAM, $(df -h $WD | awk 'NR==2{print $4}') disk free"

# --- tools (install only what is missing) ---
export PATH=$WD/software/minimap2-2.30_x64-linux:$PATH
if ! command -v bzip2 >/dev/null; then sudo apt-get -qq update && sudo apt-get -qq install -y bzip2; fi
if ! command -v minimap2 >/dev/null; then
  curl -sL https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 | tar -jxf - -C $WD/software
fi
if ! command -v samtools >/dev/null; then
  sudo apt-get -qq update && sudo apt-get -qq install -y samtools bzip2
fi
if ! command -v aws >/dev/null; then
  curl -fsSL https://awscli.amazonaws.com/v2/install.sh | sudo bash -s -- --system || sudo snap install aws-cli --classic
  hash -r
fi
echo "=== minimap2 $(minimap2 --version), $(samtools --version | head -1)"

# --- data ---
S3=s3://sg-nex-data/data
REF=reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
[ -f $REF ]     || aws s3 cp --no-sign-request --only-show-errors $S3/annotations/genome_fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa reference/
[ -f $REF.fai ] || aws s3 cp --no-sign-request --only-show-errors $S3/annotations/genome_fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai reference/
[ -f reference/Homo_sapiens.GRCh38.91.gtf ] || aws s3 cp --no-sign-request --only-show-errors $S3/annotations/gtf_file/Homo_sapiens.GRCh38.91.gtf reference/
for s in "${SAMPLES[@]}"; do
  [ -f fastq/$s.fastq.gz ] || aws s3 cp --no-sign-request --only-show-errors $S3/sequencing_data_ont/fastq/$s/$s.fastq.gz fastq/
done
ls -lh reference/ fastq/

# --- alignment ---
T=$(nproc)
MEM_GB=$(free -g | awk '/Mem:/{print $2}')

for s in "${SAMPLES[@]}"; do
  if [ -f bam/$s.bam.bai ]; then echo "=== $s already done, skipping"; continue; fi
  case $s in
    *directRNA*) flags="-ax splice -uf -k14" ;;   # direct RNA: forward strand only, smaller k-mer
    *cDNA*)      flags="-ax splice" ;;            # cDNA: either strand, default k-mer
  esac
  SPLIT=""
  if [ "$MEM_GB" -lt 28 ]; then SPLIT="-I 1g --split-prefix sam/tmp_$s"; fi

  echo "=== $(date '+%H:%M:%S')  $s"
  cmd="minimap2 $flags -t $T $SPLIT $REF fastq/$s.fastq.gz > sam/$s.sam"
  echo "\$ $cmd"
  eval "$cmd"
  echo "\$ samtools view -b sam/$s.sam > sam/$s.unsorted.bam"
  samtools view -@ $T -b sam/$s.sam > sam/$s.unsorted.bam
  echo "\$ samtools sort sam/$s.unsorted.bam -o bam/$s.bam"
  samtools sort -@ $T sam/$s.unsorted.bam -o bam/$s.bam
  echo "\$ samtools index bam/$s.bam"
  samtools index bam/$s.bam
  rm -f sam/$s.sam sam/$s.unsorted.bam
  echo "=== $(date '+%H:%M:%S')  $s done"
done
echo "ALL DONE"
ls -lh bam/

# quick summary of mapping rates (useful later for QC, question 5.3)
for s in "${SAMPLES[@]}"; do
  echo "##### $s"; samtools flagstat -@ $T bam/$s.bam | grep -E "in total|primary mapped"
done

if [ "${AUTO_STOP:-0}" = "1" ]; then echo "=== shutting down in 2 min"; sudo shutdown -h +2; fi
