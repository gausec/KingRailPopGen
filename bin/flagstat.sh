#!/bin/bash


# Mapping summary stats

for file in *.sorted.bam; do echo "$file" >> KIRAflagstats.txt; samtools flagstat -@ 30 "$file" >> KIRAflagstats.txt; echo "" >> KIRAflagstats.txt; done
echo "All done!"
