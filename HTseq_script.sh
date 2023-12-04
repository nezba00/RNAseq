#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="HTseq"
#SBATCH --time=14:00:00
#SBATCH --output=HTseq1_%j.out
#SBATCH --error=HTseq1_%j.err

# load HTseq module
module load UHTS/Analysis/HTSeq/0.9.1

sampath=/data/users/bnezar/RNA_Project/SAM_files
gtfpath=/data/users/bnezar/RNA_Project/GTF_files
home=/data/users/bnezar/RNA_Project


# Create a new directory for HTseq output
mkdir $home/HTseq

# Syntax : htseq-count [options] <alignment_files> <gff_file>
# Analyze genes
htseq-count -f sam -s reverse $sampath/1_1.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.1_HT_genes.tsv
htseq-count -f sam -s reverse $sampath/1_2.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.2_HT_genes.tsv
htseq-count -f sam -s reverse $sampath/1_3.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.3_HT_genes.tsv

htseq-count -f sam -s reverse $sampath/P_1.sam $gtfpath/all_merged.gtf > $home/HTseq/P1_HT_genes.tsv
htseq-count -f sam -s reverse $sampath/P_2.sam $gtfpath/all_merged.gtf > $home/HTseq/P2_HT_genes.tsv
htseq-count -f sam -s reverse $sampath/P_3.sam $gtfpath/all_merged.gtf > $home/HTseq/P3_HT_genes.tsv

# Analyze transcripts
htseq-count -f sam -i transcript_id -s reverse $sampath/1_1.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.1_HT_transcript.tsv
htseq-count -f sam -i transcript_id -s reverse $sampath/1_2.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.2_HT_transcript.tsv
htseq-count -f sam -i transcript_id -s reverse $sampath/1_3.sam $gtfpath/all_merged.gtf > $home/HTseq/S1.3_HT_transcript.tsv

htseq-count -f sam -i transcript_id -s reverse $sampath/P_1.sam $gtfpath/all_merged.gtf > $home/HTseq/P1_HT_transcript.tsv
htseq-count -f sam -i transcript_id -s reverse $sampath/P_2.sam $gtfpath/all_merged.gtf > $home/HTseq/P2_HT_transcript.tsv
htseq-count -f sam -i transcript_id -s reverse $sampath/P_3.sam $gtfpath/all_merged.gtf > $home/HTseq/P3_HT_transcript.tsv
