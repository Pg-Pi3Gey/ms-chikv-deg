#!/bin/bash
##  start the timer
SECONDS=0

# change working directory
cd ~/CHIKV_DEG

# STEP 1: Run fastqc #
mkdir -p ~/CHIKV_DEG/results/fastqc/untrimmed
cd ~/CHIKV_DEG/raw_seq/untrimmed_seq
mkdir -p ~/CHIKV_DEG/results/fastqc/untrimmed
fastqc -t 16 *.fastq.gz -o ~/CHIKV_DEG/results/fastqc/untrimmed
echo "|> Quality Check via FASTQC on untrimmed sequencence completed <|"
# create single report via multiqc
multiqc ~/CHIKV_DEG/results/fastqc/untrimmed -o ~/CHIKV_DEG/results/multiqc -n untrimmed_multiqc_report.html

##  Total runtime
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."