#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="five prime intersect"
#SBATCH --time=00:30:00
#SBATCH --output=./output_files/five_intersect_%j.out
#SBATCH --error=./error_files/five_intersect%j.err


# Load bedtools
module load UHTS/Analysis/BEDTools/2.29.2

home=/data/users/bnezar/RNA_Project/step_6
ref=/data/courses/rnaseq_course/lncRNAs/Project1/references
#-->refTSS_v4.1_human_coordinate.hg38.bed

cd $home

# Use bedtools to find regions which do intersect with TSS
# Options:
#   s: only report hits on the same strand
#   wa: write entry from a
#   wb: write entry from b
#   a: Each feature in A is compared to B in search of overlaps (BAM/BED/GFF/VCF)
#   b: One or more BAM/BED/GFF/VCF file(s) “B”
bedtools intersect -s -wa -wb -a five_prime.bed \
                              -b $ref/refTSS_v4.1_human_coordinate.hg38.bed \
                               > $home/output_files/tss_intersect.bed


