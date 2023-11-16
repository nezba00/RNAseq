#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G

# Load FastQC
module load UHTS/Quality_control/fastqc/0.11.9;

# Count how many reads there are in each parental file
FASTQdir="/data/courses/rnaseq_course/lncRNAs/fastq"

#for file in "$FASTQdir"/P*.gz;
#do 
#    lines=$(zcat "$file" | wc -l)
#    name=$(basename "$file")
#    # Check if $lines is bigger than 0
#    if [ "$lines" -gt 0 ]; then
#        result=$(($lines / 4))
#        echo "$name: $result" >> /data/users/bnezar/Amount-Replicates.txt
#    else 
#        echo "$name: Lines == 0" >> /data/users/bnezar/Amount-Replicates.txt
#    fi
#done


# Generate fastqcs of all sequences
for file in "$FASTQdir"/P*.gz;
do
    fastqc -o "/data/users/bnezar/RNA_Project/FastQC_Files" "$file"
done

