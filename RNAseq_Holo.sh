#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G

# Holoclonal
# Clones 1.1/1.2/1.5

# Home Directory
# /data/users/bnezar

#wc -l /data/courses/rnaseq_course/lncRNAs/fastq
HOME=/data/users/bnezar
FASTQdir="/data/courses/rnaseq_course/lncRNAs/fastq"

# Count how many reads there are in each holo file
for file in "$FASTQdir"/1*.gz;
do 
    lines=$(zcat "$file" | wc -l)
    name=$(basename "$file")
    # Check if $lines is bigger than 0
    if [ "$lines" -gt 0 ]; then
        result=$(($lines / 4))
        echo "$name: $result" >> /data/users/bnezar/Amount-Replicates.txt
    else 
        echo "$name: Lines == 0" >> /data/users/bnezar/Amount-Replicates.txt
    fi
done

# Count how many reads there are in each parental file
for file in "$FASTQdir"/P*.gz;
do 
    lines=$(zcat "$file" | wc -l)
    name=$(basename "$file")
    # Check if $lines is bigger than 0
    if [ "$lines" -gt 0 ]; then
        result=$(($lines / 4))
        echo "$name: $result" >> /data/users/bnezar/Amount-Replicates.txt
    else 
        echo "$name: Lines == 0" >> /data/users/bnezar/Amount-Replicates.txt
    fi
done


# Generate fastqcs of all sequences
for file in "$FASTQdir"/1*.gz;
do
    fastqc -o "$HOME/RNA_Project/FastQC_Files" "$file"
done
