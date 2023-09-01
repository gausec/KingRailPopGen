#!/bin/bash

# Check for the input file argument
if [ $# -eq 0 ]; then
  echo "Usage: $0 <input_file.thetas.idx.pestPG>"
  exit 1
fi

# Input file name
input_file="$1"

# Output file names
tW_file="tW.txt"
cleaned_tW_file="cleaned_tW.txt"
average_tW_file="Average_tW.txt"

# Calculate Watterson's theta (tW)
awk '{print $4 / $14}' "$input_file" >> "$tW_file"

# Remove empty rows
awk '$1 != "-nan"' "$tW_file" > "$cleaned_tW_file"

# Calculate Waterson's theta by taking the average
awk '{ sum += $1 } END { print sum / NR }' "$cleaned_tW_file" > "$average_tW_file"

# Print the input file name and final value to the output file
echo "Input File: $input_file" >> "$average_tW_file"

# Display completion message
echo "Process is complete. Results in Average_tW.txt"

