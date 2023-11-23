#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name="Stringtie"
#SBATCH --time=10:00:00

# Load string tie module
module load UHTS/Aligner/stringtie/1.3.3b

home_path=/data/users/bnezar/RNA_Project/
reference_path=/data/courses/rnaseq_course/lncRNAs/Project1/references

# create a directory for gtf files
mkdir $home_path/GTF_files

# stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
# Template guided transcriptome assembly for all samples
# -p to specify amount of cpus
stringtie -o $home_path/GTF_files/1_1_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_1_1.bam
stringtie -o $home_path/GTF_files/1_2_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_1_2.bam
stringtie -o $home_path/GTF_files/1_3_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_1_3.bam

stringtie -o $home_path/GTF_files/P_1_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_P_1.bam
stringtie -o $home_path/GTF_files/P_2_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_P_2.bam
stringtie -o $home_path/GTF_files/P_3_transcripts -p 8 --rf -G $reference_path/gencode.v44.chr_patch_hapl_scaff.annotation.gtf  $home_path/BAM_files/sorted_P_3.bam



