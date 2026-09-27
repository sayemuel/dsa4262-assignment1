#!/bin/bash
# Overnight chain: Task 4 Bambu -> Task 5 scenario 1 -> scenario 2 (-resume) -> shutdown
source ~/miniforge3/etc/profile.d/conda.sh && conda activate genomics
sudo shutdown -h +480 "backstop"          # hard stop after 8 h no matter what

# Task 4 (skip if a finished result already exists)
cd ~/task4
if ! grep -q "^DONE" bambu.log 2>/dev/null; then
  Rscript run_bambu.R > bambu.log 2>&1
fi

# Task 5 only if the small test passed
if grep -q TEST_PASSED ~/task5_test/test.log 2>/dev/null; then
  cd ~/task5
  REF=~/task3/reference
  S1="nextflow run main.nf -with-report report_scenario1.html -with-timeline timeline_scenario1.html --reads '$HOME/task3/fastq/*.fastq.gz' --refFa $REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --refGtf $REF/Homo_sapiens.GRCh38.91.gtf --outdir results_scenario1"
  S2="nextflow run main.nf -resume -with-report report_scenario2.html -with-timeline timeline_scenario2.html --reads '$HOME/task3/fastq/*.fastq.gz' --refFa $REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --refGtf $REF/Homo_sapiens.GRCh38.91.gtf --outdir results_scenario2 --annotations false"
  # 'script' records the real terminal output so it can be replayed tomorrow with: cat scenario1.typescript
  script -q -c "echo '\$ $S1'; $S1" scenario1.typescript
  script -q -c "echo '\$ $S2'; $S2" scenario2.typescript
fi
echo "OVERNIGHT DONE $(date)" >> ~/overnight.log
sudo shutdown -h +2
