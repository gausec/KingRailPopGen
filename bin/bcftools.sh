#!/bin/bash

# variant calling pipeline
bcftools mpileup --threads 40 -f Rallus_crepitans_1.0.fasta -Ou -q 20 -Q 20 -C 50 -I -d 1000000 -a AD 43617.sorted.bam | bcftools call --threads 40 -mv -Ov -f GQ -o 43617_output_raw.vcf

echo "done"

