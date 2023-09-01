#!/bin/bash

# set working directory
bam_directory="/home/gausec21/WholeGenomes/SortedBAM"

# loop through *.sorted.bam files
for bam_file in "$bam_directory"/*.sorted.bam; do
    # Extract the filename without the extension
    file_prefix=$(basename "$bam_file" .sorted.bam)

    # bcftools commands
    bcftools mpileup -Ou -f CLRAindex/Rallus_crepitans_1.0.fasta --threads 20 "$bam_file" | bcftools call -mv -Ob -o "${file_prefix}_output_filtered.bcf"

done
