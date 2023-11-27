#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="Kall_Index"
#SBATCH --time=05:00:00
#SBATCH --output=kallIndex_%j.out
#SBATCH --error=kallIndex_%j.err

# Load Kallisto module
module load UHTS/Analysis/kallisto/0.46.0

# Define path variables
home=/data/users/bnezar/RNA_Project
reference=/data/courses/rnaseq_course/lncRNAs/Project1/references

#Create Kallisto directory
mkdir kallisto

# Run kallisto index on reference file
kallisto index --index="$home/kallisto/GRCh38_kallisto_index" $reference/GRCh38.genome.fa

