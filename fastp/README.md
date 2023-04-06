# Preprocessing for FastQ files
This repository contains steps for the preprocessing of paired-end whole genome data of 4 species (king rail, clapper rail, Virginia rail, and eastern bluebird). All steps below occured in the command line once the original WGS data was saved to an external hard drive. 
---
### Step 1: File transfer, extraction, and organization
1. Create a subdirectory within the remote server
```
mkdir WholeGenomes
```
2. Transfer tar files to remote server using the `scp` command
3. Extract folder containing zipped fastq files 
```
tar -xvf sequences.tar
```
4. Navigate to the extracted directory using the `cd` command, and check that the fastq.gz files were successfully extracted
5. Navigate back to previous directory using `cd ..` command and remove the tar files to conserve server space
```
rm sequences.tar
```
6. Nvaigate to forward reads folder to move all forward and reverse sequences to the parent directory for easier processing. Repeat for reverse reads. Remove empty subdirectores. 
```
 mv *.gz ..
```
7. Create a subdirectory for the cleaned reads
```
mkdir cleaned
```
*Note: The owner of the Biology Department server installed fastp for this project.*

### Step 2: Quality control
1. Generate lists of all forward and reverse reads for peace of mind:
```
for FileName in *_R1_001.fastq.gz; do ls $FileName; done 
for FileName in *_R2_001.fastq.gz; do ls $FileName; done 
```
2. Use the `cut` command to extract the first field of the filename before the first _ character, store it as the variable `name`, and print it out for each file using the `echo` command. I'm doing this because I want to rename my output files by their shorter sample names.
```
for FileName in *_R1_001.fastq.gz; do name=`ls $FileName | cut -d"_" -f1`; echo $name; done 
```
3. Use a loop to run the `fastp` command (without changing default quality parameters) for all fastq files. This command specifies the input files for both the forward (`-i`) and reverse (`-I`) reads, as well as the output files for the processed reads (`-o` and `-O`). The -`-html` and `--json` options specify the output files for the HTML and JSON reports, respectively. *Note: Fastp is a tool that can be used for quality control and trimming of high-throughput sequencing data in fastq format, including fastq.gz compressed files: https://github.com/OpenGene/fastp* 
```
for FileName in *_R1_001.fastq.gz; do name=`ls $FileName | cut -d"_" -f1`; fastp -i $name*_R1_001.fastq.gz -I $name*_R2_001.fastq.gz -o cleaned/$FileName-R1.fastq.gz -O cleaned/$FileName-R2.fastq.gz --html $FileName.html --json $FileName.json ; done 
```


