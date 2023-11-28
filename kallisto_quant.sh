#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="Kall_Quant"
#SBATCH --time=05:00:00
#SBATCH --output=kallQuant_%j.out
#SBATCH --error=kallQuant_%j.err

# load kallisto module
module load UHTS/Analysis/kallisto/0.46.0

home=/data/users/bnezar/RNA_Project
reference=/data/courses/rnaseq_course/lncRNAs/Project1/references

# Run Kallisto quant
# Args:
# -i: index file
# -o: output directory
# --rf-stranded: for fr-firststrand libraries (TruSeq stranded sample prep kits)
# -t: Number of threads to use
# Usage: kallisto quant [arguments] FASTQ-files
kallisto quant -t 8 --rf-stranded \
                -i $home/kallisto/GRCh38_kallisto_index \
                -o $home/kallisto \
                "1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz" "1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz" \
                "1_2_L3_R1_001_DnNWKUYhfc9S.fastq.gz" "1_2_L3_R2_001_SNLaVsTQ6pwl.fastq.gz \
                "1_5_L3_R1_001_iXvvRzwmFxF3.fastq.gz" "1_5_L3_R2_001_iXCMrktKyEh0.fastq.gz \
                "P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz" "P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz" \
                "P2_L3_R1_001_R82RphLQ2938.fastq.gz" "P2_L3_R2_001_06FRMIIGwpH6.fastq.gz" \
                "P3_L3_R1_001_fjv6hlbFgCST.fastq.gz" "P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz"


