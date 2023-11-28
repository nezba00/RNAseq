#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="Kall_Index"
#SBATCH --time=05:00:00
#SBATCH --output=kallIndex_%j.out
#SBATCH --error=kallIndex_%j.err

# Load Kallisto and cufflinks module
module load UHTS/Analysis/kallisto/0.46.0
module load UHTS/Assembler/cufflinks/2.2.1;

# Define path variables
home=/data/users/bnezar/RNA_Project
reference=/data/courses/rnaseq_course/lncRNAs/Project1/references

#Create Kallisto directory
mkdir $home/kallisto

# Create a fasta with all transcripts
gffread -g $reference/GRCh38.genome.fa -w $home/kallisto/all_transcripts.fa $home/GTF_files/all_merged.gtf

# Run kallisto index on transcripts fasta file
kallisto index -i $home/kallisto/GRCh38_kallisto_index $home/kallisto/all_transcripts.fa

