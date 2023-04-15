#!/bin/bash

echo "Creating index"

# Create index directory & index files from the clapper rail genome
mkdir index
mv *.fasta index/
bwa index index/*.fasta

echo "Aligning"

# Assign basename & align the trimmed reads to the reference genome using the read alignment tool
mkdir SAM
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 16 index/Rallus_crepitans_1.0.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done &

echo "SAM > BAM"
# Convert SAM to BAM using samtools (this saves a lot of space by converting to binary data)
for r1file in *.sam; do base=$(basename $r1file .sam); samtools view -@ 40 $r1file > ${base}.bam; done

echo "Sorting"

# Arrange the aligned reads in the BAM files based on their coordinates in the reference genome (This is called sorting)
for file in *.bam; do echo "Sorting ${file}..."; samtools sort -@ 40 -T temp_dir ${file} -o sorted/${file%.*}.sorted.bam; done*

echo "All done!"

