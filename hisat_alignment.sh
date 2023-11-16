#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00

# Load Hisat2 module
module load UHTS/Aligner/hisat/2.2.1

home_path=/data/users/bnezar/RNA_Project
references_path=/data/courses/rnaseq_course/lncRNAs/fastq

# syntax for hisat1: hisat2 -x "reference" {-1 data_1 -2 data_2}

# Alignment for test strains
# --phred33 NOt sure if default but better to specify
# -p 8 to make use of all 8 cores
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz  -2 $references_path/1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz -S $home_path/SAM/1_1.sam
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/1_2_L3_R1_001_DnNWKUYhfc9S.fastq.gz -2 $references_path/1_2_L3_R2_001_SNLaVsTQ6pwl.fastq.gz -S $home_path/SAM/1_2.sam
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/1_5_L3_R1_001_iXvvRzwmFxF3.fastq.gz  -2 $references_path/1_5_L3_R2_001_iXCMrktKyEh0.fastq.gz -S $home_path/SAM/1_3.sam

# Alignment for parental strains
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz  -2 $references_path/P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz -S $home_path/SAM/P_1.sam
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/P2_L3_R1_001_R82RphLQ2938.fastq.gz  -2 $references_path/P2_L3_R2_001_06FRMIIGwpH6.fastq.gz -S $home_path/SAM/P_2.sam
hisat2 --phred33 -p 8 -x $home_path/Reference_index -1 $references_path/P3_L3_R1_001_fjv6hlbFgCST.fastq.gz  -2 $references_path/P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz -S $home_path/SAM/P_3.sam
