#!/bin/bash
###  start the timer
SECONDS=0

# change working directory
cd ~/CHIKV_DEG


## STEP 1: Qualitycheck via FASTQC ##
mkdir -p ~/CHIKV_DEG/results/fastqc/untrimmed
fastqc -t 14 ~/CHIKV_DEG/raw_seq/untrimmed_seq/*.fastq.gz -o ~/CHIKV_DEG/results/fastqc/untrimmed
echo "|> Quality Check via FASTQC on untrimmed sequencence completed <|"
# create single report via multiqc #
multiqc ~/CHIKV_DEG/results/fastqc/untrimmed -o ~/CHIKV_DEG/results/multiqc -n untrimmed_multiqc_report.html --interactive


## STEP 2: Trimming bad reads via Trimmomatic ##
cd ~/CHIKV_DEG/raw_seq/untrimmed_seq/
cp ~/CHIKV_DEG/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa .
for infile in *_1.fastq.gz
do
    base=$(basename ${infile} _1.fastq.gz)
    java -jar ~/CHIKV_DEG/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 14 -phred33 ${infile} ${base}_2.fastq.gz \
    ${base}_1_paired.fastq.gz ${base}_1_unpaired.fastq.gz \
    ${base}_2_paired.fastq.gz ${base}_2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True LEADING:7 TRAILING:7 SLIDINGWINDOW:5:20 MINLEN:35
done
# move trimmed sequence to different folder #
mkdir ../trimmed_seq
mkdir ../unpaired_seq
mv *unpaired* ../unpaired_seq
mv *pair* ../trimmed_seq
echo "|> Trimming completed! <|"

# Quality check the trimmed sequences via FASTQC #
mkdir -p ~/CHIKV_DEG/results/fastqc/trimmed
fastqc -t 14 ~/CHIKV_DEG/raw_seq/trimmed_seq/*.fastq.gz -o ~/CHIKV_DEG/results/fastqc/trimmed
# create single report via multiqc #
multiqc --interactive ~/CHIKV_DEG/results/fastqc/trimmed -o ~/CHIKV_DEG/results/multiqc -n trimmed_multiqc_report.html

# delete raw & unpaired reads to save storage #
rm ~/CHIKV_DEG/raw_seq/trimmed_seq/*_unpaired*
rm -rf ~/CHIKV_DEG/raw_seq/untrimmed_seq

## STEP 3: Alignments via HISAT2 ##

# # get the genome index #
# mkdir -p ~/CHIKV_DEG/tools/hisat2
# cd ~/CHIKV_DEG/tools/hisat2
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# tar -xvzf grch38_genome.tar.gz && rm grch38_genome.tar.gz

# align the trimmed reads to the genome #
cd ~/CHIKV_DEG/raw_seq/trimmed_seq/
mkdir -p ~/CHIKV_DEG/results/hisat2
for i in `ls *.fastq.gz | sed 's/_[12]_paired.fastq.gz//g' | sort -u`
do
echo "Aligning: $i" && hisat2 -p 14 --rna-strandness RF -x ~/CHIKV_DEG/tools/hisat2/grch38/genome -1 ${i}_1_paired.fastq.gz -2 ${i}_2_paired.fastq.gz | samtools sort -o ~/CHIKV_DEG/results/hisat2/${i}.hisat.sorted.bam
done
echo "|> Alignment Completed! <|"

## STEP 4: Quantification via Featurecount ##
mkdir -p ~/CHIKV_DEG/tools/featurecount
cd ~/CHIKV_DEG/tools/featurecount

# # # get the gtf file #
# # wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
# # gzip -d Homo_sapiens.GRCh38.109.gtf.gz

# run featurecounts #
mkdir -p ~/CHIKV_DEG/results/featurecounts/
featureCounts -T 14 -s 2 -p --countReadPairs -a Homo_sapiens.GRCh38.109.gtf -o ~/CHIKV_DEG/results/featurecounts/featurecounts_output.txt ../../results/hisat2/*.bam
echo "|> Quantification Completed! <|"

# create single report via multiqc #
multiqc --interactive ~/CHIKV_DEG/results/. -o ~/CHIKV_DEG/results/multiqc -n overall_multiqc_report.html

###  Total runtime
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."