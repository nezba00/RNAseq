#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00

# Load Trimmomatic module
module load Trimmomatic/0.39-Java-11;

home_path=/data/users/bnezar/RNA_Project
references_path=/data/courses/rnaseq_course/lncRNAs/fastq

trimmomatic PE -threads 8 $references_path/1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz \
                          $references_path/1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz \
                          $home_path/trimmomatic_out/1_1_R1.trimmed.fastq \
                          $home_path/trimmomatic_out/1_1_R1un.trimmed.fastq \
                          $home_path/trimmomatic_out/1_1_R2.trimmed.fastq \
                          $home_path/trimmomatic_out/1_1_R2un.trimmed.fastq \
                          ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20

java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36