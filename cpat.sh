#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="CPAT"
#SBATCH --time=02:00:00
#SBATCH --output=./output_files/CPAT_%j.out
#SBATCH --error=./error_files/CPAT_%j.err

# Load CPAT and bedtools
module load SequenceAnalysis/GenePrediction/cpat/1.2.4
module load UHTS/Analysis/BEDTools/2.29.2

home=/data/users/bnezar/RNA_Project/step_6
ref=/data/courses/rnaseq_course/lncRNAs/Project1/references

cd $home

# Get fasta references to bed file with novel transcripts
#   fi: fasta in
#   bed: bed reference
#   fo: specify output dir
bedtools getfasta -name+ -s -fi $ref/GRCh38.genome.fa -bed $home/novel.bed \
                    -fo $home/output_files/novel.fa

## Run CPAT to predict coding potential
# Human_Hexamer.tsv from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/
# logit model from same place
# OPTIONS:
#   x: The hexamer frequency table. The prebuilt tables for Human, Mouse, Fly, Zebrafish are available
#   --antisense: Also search for ORFs from the anti-sense strand. 
#   d: Logistic regression model. The prebuilt models for Human, Mouse, Fly, Zebrafish are available
#   --top-orf: Number of ORF candidates --> this is default
#   g: Genomic sequnence(s) of RNA in FASTA or BED
#   r: Reference genome sequences in FASTA format --> important when g is a .bed
#   o: Output file
cpat.py -x $ref/Human_Hexamer.tsv -d $ref/Human_logitModel.RData \
        -g $home/output_files/novel.fa \
        -o $home/output_files/cpat_out



