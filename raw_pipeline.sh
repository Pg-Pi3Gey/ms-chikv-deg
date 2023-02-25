#!/bin/bash
###  start the timer
SECONDS=0

# change working directory
cd ~/CHIKV_DEG

## STEP 1: Qualitycheck via FASTQC ##
mkdir -p ~/CHIKV_DEG/results/fastqc/untrimmed
fastqc -t 16 ~/CHIKV_DEG/raw_seq/untrimmed_seq/*.fastq.gz -o ~/CHIKV_DEG/results/fastqc/untrimmed
echo "|> Quality Check via FASTQC on untrimmed sequencence completed <|"
# create single report via multiqc #
multiqc ~/CHIKV_DEG/results/fastqc/untrimmed -o ~/CHIKV_DEG/results/multiqc -n untrimmed_multiqc_report.html


## STEP 2: Trimming bad reads via Trimmomatic ##
cd ~/CHIKV_DEG/raw_seq/untrimmed_seq/
cp ~/CHIKV_DEG/tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa .
for infile in *_1.fastq.gz
do
    base=$(basename ${infile} _1.fastq.gz)
    java -jar ~/CHIKV_DEG/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ${infile} ${base}_2.fastq.gz \
    ${base}_1_paired.fastq.gz ${base}_1_unpaired.fastq.gz \
    ${base}_2_paired.fastq.gz ${base}_2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 MINLEN:35
done
# move trimmed sequence to different folder #
mkdir ../trimmed_seq
mv *pair* ../trimmed_seq
echo "|> Trimmomatic finished running! <|"

# Quality check the trimmed sequences via FASTQC #
mkdir -p ~/CHIKV_DEG/results/fastqc/trimmed
fastqc -t 16 ~/CHIKV_DEG/raw_seq/trimmed_seq/*.fastq.gz -o ~/CHIKV_DEG/results/fastqc/trimmed
# create single report via multiqc #
multiqc ~/CHIKV_DEG/results/fastqc/trimmed -o ~/CHIKV_DEG/results/multiqc -n trimmed_multiqc_report.html


# # delete unpaired reads to save storage #
# rm ~/CHIKV_DEG/raw_seq/trimmed_seq/*_unpaired*
# rm -rf ~/CHIKV_DEG/raw_seq/untrimmed_seq

##  Total runtime
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."