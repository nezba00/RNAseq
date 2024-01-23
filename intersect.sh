#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="bdtools intersect"
#SBATCH --time=00:30:00
#SBATCH --output=./output_files/intersect_%j.out
#SBATCH --error=./error_files/intersect%j.err

# filter bed files for transcripts in intergenic regions

# Load bedtools
module load UHTS/Analysis/BEDTools/2.29.2

home=/data/users/bnezar/RNA_Project/step_6

cd $home
# Use bedtools to find regions which do not intersect with known genes
# Options:
#   v: Only report those entries in A that have no overlap in B
#   a: Each feature in A is compared to B in search of overlaps (BAM/BED/GFF/VCF)
#   b: One or more BAM/BED/GFF/VCF file(s) “B”
bedtools intersect -v -a novel.bed -b known.bed > $home/output_files/intersect.bed