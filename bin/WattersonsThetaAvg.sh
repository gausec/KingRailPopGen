#!/bin/bash

# Note for future me: = and == are for string comparisons while -eq is for numeric comparisons
# The "$#" stores the total number of arguments

# check for input file
if [ $# -eq 0 ]; then
  echo "Usage: $0 <input_file.thetas.idx.pestPG>"
  exit 1
fi

# input file name
input_file="$1"

# output file names
tW_file="tW.txt"
cleaned_tW_file="cleaned_tW.txt"
average_tW_file="Average_tW.txt"

# calculate Watterson's theta 
awk '{print $4 / $14}' "$input_file" >> "$tW_file"

# remove empty rows
awk '$1 != "-nan"' "$tW_file" > "$cleaned_tW_file"

# calculate avg. Waterson's theta
awk '{ sum += $1 } END { print sum / NR }' "$cleaned_tW_file" > "$average_tW_file"

# print the input file name + final value to output file
echo "Input File: $input_file" >> "$average_tW_file"

# done
echo "Process is complete. Results in Average_tW.txt"

