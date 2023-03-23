# Whole_Genome_QC
This repository contains files and information pertaining to quality control for whole genome data of 4 species (king rail, clapper rail, Virginia rail, and eastern bluebird). All steps below occured once the original WGS data was saved to an external hard drive. 
---
### Step 1: File transfer & extraction
1. Create a subdirectory within the remote server
  * `mkdir WholeGenomes`
2. Transfer tar files to remote server using the `scp` command
3. Extract directory containing zipped fastq files 
  * `tar -xvf sequences.tar`
4. Navigate to the extracted directory using the `cd` command, and check that the fastq.gz files were successfully extracted
5. Navigate back to previous directory using `cd ..` command and remove the tar files to conserve server space
  * `rm sequences.tar`

### Step 2: Quality control
1. In the subdirectory containing the fastq.gz files, use the `fastp` command
