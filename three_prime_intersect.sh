#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="three prime intersect"
#SBATCH --time=00:30:00
#SBATCH --output=./output_files/three_intersect_%j.out
#SBATCH --error=./error_files/three_intersect%j.err

# This works analog to five_prime_intersect

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
bedtools intersect -s -wa -wb -a three_prime_light.bed \
                              -b $ref/atlas.clusters.2.0.GRCh38.96.bed \
                               > $home/output_files/three_prime_intersect.bed