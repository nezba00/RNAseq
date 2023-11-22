#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --job-name="SAM_BAM Convert"
#SBATCH --time=20:00:00


# load SAMtools
module load UHTS/Analysis/samtools/1.10

# Create directory for BAM files
mkdir /data/users/bnezar/RNA_Project/BAM

# convert .sam to .bam for all samples
# @ to specify cores
# -Sb to specify sam to bam
# -o to specify output file
samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/1_1.bam /data/users/bnezar/RNA_Project/SAM/1_1.sam
samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/1_2.bam /data/users/bnezar/RNA_Project/SAM/1_2.sam
samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/1_3.bam /data/users/bnezar/RNA_Project/SAM/1_3.sam

samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/P_1.bam /data/users/bnezar/RNA_Project/SAM/P_1.sam
samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/P_2.bam /data/users/bnezar/RNA_Project/SAM/P_2.sam
samtools view -@ 8 -Sb -o /data/users/bnezar/RNA_Project/BAM/P_3.bam /data/users/bnezar/RNA_Project/SAM/P_3.sam

# Sort the BAM files
samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_1_1.bam /data/users/bnezar/RNA_Project/BAM/1_1.bam
samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_1_2.bam /data/users/bnezar/RNA_Project/BAM/1_2.bam
samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_1_3.bam /data/users/bnezar/RNA_Project/BAM/1_3.bam

samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_P_1.bam /data/users/bnezar/RNA_Project/BAM/P_1.bam
samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_P_2.bam /data/users/bnezar/RNA_Project/BAM/P_2.bam
samtools sort -O bam -o /data/users/bnezar/RNA_Project/BAM/sorted_P_3.bam /data/users/bnezar/RNA_Project/BAM/P_3.bam


# Index the files
samtools index /data/users/bnezar/RNA_Project/BAM/sorted_1_1.bam
samtools index /data/users/bnezar/RNA_Project/BAM/sorted_1_2.bam
samtools index /data/users/bnezar/RNA_Project/BAM/sorted_1_3.bam

samtools index /data/users/bnezar/RNA_Project/BAM/sorted_P_1.bam
samtools index /data/users/bnezar/RNA_Project/BAM/sorted_P_2.bam
samtools index /data/users/bnezar/RNA_Project/BAM/sorted_P_3.bam