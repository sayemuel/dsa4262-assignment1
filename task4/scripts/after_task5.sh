#!/bin/bash
# Run Task 4 Bambu after the Task 5 chain finishes, then shut down
until grep -q "OVERNIGHT DONE" ~/overnight.log 2>/dev/null; do sleep 20; done
sudo shutdown -c                      # cancel the chain's "shutdown in 2 min"
sudo shutdown -h +300 "backstop"      # new backstop: 5 h
source ~/miniforge3/etc/profile.d/conda.sh && conda activate genomics
cd ~/task4
NCORE=4 Rscript run_bambu.R > bambu.log 2>&1
grep -q "^DONE" bambu.log || NCORE=1 Rscript run_bambu.R > bambu.log 2>&1   # fallback: one sample at a time
echo "TASK4 FINISHED $(date)" >> ~/overnight.log
sudo shutdown -h +2
