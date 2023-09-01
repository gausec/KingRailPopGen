#!/bin/bash


# Note for future me: = and == are for string comparisons while -eq is for numeric comparisons
# The "$#" stores the total number of arguments


# Check for the input file
if [ $# -eq 0 ]; then
  echo "Usage: $0 <input_file.thetas.idx.pestPG>"
  exit 1
fi

# Input file name
input_file="$1"

# output file names
pi_file="pi.txt"
cleaned_pi_file="cleaned_pi.txt"
average_pi_file="Average_pi.txt"

# Extract theta P (tP) col & divide by the number of sites (nSites row)
awk '{print $5 / $14}' "$input_file" >> "$pi_file"

# remove empty rows
awk '$1 != "-nan"' "$pi_file" > "$cleaned_pi_file"

# Calculate pi by taking the average
awk '{ sum += $1 } END { print sum / NR }' "$cleaned_pi_file" > "$average_pi_file"

# Print input file name & final value to output file
echo "Input File: $input_file" >> "$average_pi_file"

# End
echo "Calculation complete. Results in Average_pi.txt"
