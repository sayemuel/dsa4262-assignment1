#!/bin/bash
# Wait for Task 4 Bambu, then rerun Task 5 scenario 2 (no annotations) with -resume
until grep -qE "^DONE|Execution halted" ~/task4/bambu.log; do sleep 30; done
source ~/miniforge3/etc/profile.d/conda.sh && conda activate genomics
cd ~/task5
rm -rf results_scenario2 report_scenario2.html timeline_scenario2.html scenario2.typescript
REF=~/task3/reference
S2="nextflow run main.nf -resume -with-report report_scenario2.html -with-timeline timeline_scenario2.html --reads '$HOME/task3/fastq/*.fastq.gz' --refFa $REF/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --refGtf $REF/Homo_sapiens.GRCh38.91.gtf --outdir results_scenario2 --annotations false"
script -q -c "echo '\$ $S2'; $S2" scenario2.typescript
echo "SCENARIO2 RERUN FINISHED" > ~/task5/s2_done.flag
