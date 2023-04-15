#!/bin/bash

# This is a bash for loop that iterates over all .bam files in the current directory & calculates the average depth of coverage using samtools depth and awk

for file in *.bam; do echo "$file" >> AvgDepth.txt; samtools depth $file | awk '{sum+=$3} END {print "Average depth:", sum/NR}' > "${file%.bam}_depth.txt";  echo "" >> AvgDepth.txt; done


 
