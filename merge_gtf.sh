#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="merge"
#SBATCH --time=10:00:00


# Load string tie module
module load UHTS/Aligner/stringtie/1.3.3b

home_path=/data/users/bnezar/RNA_Project
reference_path=/data/courses/rnaseq_course/lncRNAs/Project1/references


# Concatenate all GTF file names into a merged .txt file
ls $home_path/GTF_files/* >> $home_path/GTF_files/all_merged.txt

# Generate merged GTF
stringtie --merge -p 8 --rf -o $home_path/GTF_files/all_merged.gtf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf $home_path/GTF_files/all_merged.txt




