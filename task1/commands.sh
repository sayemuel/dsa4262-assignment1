# Task 1 - Unix tutorial
# run in the folder containing Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

# a)
ls -lah

# b)
cd
cd ~
cd $HOME

# c) 194 sequences (25 chromosomes + 169 scaffolds)
grep -c ">" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
grep -c "dna_sm:chromosome" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

# d)
head -n 500 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa > chr1_500_lines.fa

# e)
grep -A 500 ">" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa > GRCH38_500_lines.fa

# f) ~4.6 MB
grep -m 75000 "ACCACCT" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa > sequence_hits.txt
ls -lh sequence_hits.txt

# g) ~1.4 MB
gzip sequence_hits.txt
ls -lh sequence_hits.txt.gz

# h)
grep -m 75000 "ACCACCT" Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa | gzip > sequence_hits.txt.gz
