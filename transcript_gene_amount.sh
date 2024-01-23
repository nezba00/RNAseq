#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="Counts"
#SBATCH --time=00:05:00
#SBATCH --output=./output_files/counts_%j.out
#SBATCH --error=./error_files/counts_%j.err

kallisto_dir=/data/users/bnezar/RNA_Project/kallisto
kall_output=/data/users/bnezar/RNA_Project/kallisto/sleuth_in

cd $kall_output

# Create Summary file
touch $kallisto_dir/summary.txt

for replicate in * ; do
    cd $kall_output/$replicate

    # calculate sum of tpf
    # -F for separator
    tpf_sum=$(awk -F'\t' 'NR > 1 {sum += ($5)} END {print sum}' abundance.tsv)

    # count transcripts
    # -f for column or column range
    transcripts=$(awk -F"\t" 'NR>1 && $4 > 0 {++count} END{print count}' abundance.tsv)

    # novel transcripts
    # same as transcripts but we only want MSTR and not enst
    novel_transcripts=$(awk -F"\t" 'NR>1 && $1 ~ /^MSTR/ && $4 > 0 {++count} END{print count}' abundance.tsv)

    # count genes
    # same as transcripts but the first field should end with 1
    genes=$(awk -F"\t" 'NR>1 && $1 ~ /1$/ && $4 > 0 {++count} END{print count}' abundance.tsv)

    # count novel genes
    # same as for genes but first row should start with MSTR
    novel_genes=$(awk -F"\t" 'NR>1 && $1 ~ /^MSTR.*1$/ && $4 > 0 {++count} END{print count}' abundance.tsv)

    # Write information to Summary file
    echo "Replicate: $replicate" >> $kallisto_dir/summary.txt
    echo "TPF sum: $tpf_sum" >> $kallisto_dir/summary.txt
    echo "Number of Transcripts: $transcripts" >> $kallisto_dir/summary.txt
    echo "Number of novel Transcripts: $novel_transcripts" >> $kallisto_dir/summary.txt
    echo "Number of genes: $genes" >> $kallisto_dir/summary.txt
    echo "Number of novel genes: $novel_genes" >> $kallisto_dir/summary.txt
done    











