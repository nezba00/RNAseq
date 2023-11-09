#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00

home_path=/data/users/bnezar/RNA_Project
references_path=/data/courses/rnaseq_course/lncRNAs/Project1/references


module load UHTS/Aligner/hisat/2.2.1

hisat2-build -p8 -f $references_path/GRCh38.genome.fa $home_path/Reference_index/GRCh38


