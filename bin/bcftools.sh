#!/bin/bash

# Set directory
bam_directory="/home/gausec21/WholeGenomes/SortedBAM"

# Loop through each *.sorted.bam file
for bam_file in "$bam_directory"/*.sorted.bam; do
    # Extract the filename without the extension
    file_prefix=$(basename "$bam_file" .sorted.bam)

    # Perform the bcftools commands using the current bam_file and file_prefix
    bcftools mpileup -Ou -f CLRAindex/Rallus_crepitans_1.0.fasta --threads 20 "$bam_file" | bcftools call -mv -Ob -o "${file_prefix}_output_filtered.bcf"

done
