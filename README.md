# Whole_Genome_QC
This repository contains files and information pertaining to quality control for paired-end whole genome data of 4 species (king rail, clapper rail, Virginia rail, and eastern bluebird). All steps below occured in the command line once the original WGS data was saved to an external hard drive. 
---
### Step 1: File transfer & extraction
1. Create a subdirectory within the remote server
```
mkdir WholeGenomes
```
2. Transfer tar files to remote server using the `scp` command
3. Extract directory containing zipped fastq files 
```
tar -xvf sequences.tar
```
4. Navigate to the extracted directory using the `cd` command, and check that the fastq.gz files were successfully extracted
5. Navigate back to previous directory using `cd ..` command and remove the tar files to conserve server space
```
rm sequences.tar
```

### Step 2: Quality control
1. First, since my forward and reverse reads were provided in separate directories, I will generate a list of all reads:
```
ls forward_reads/*.fq.gz > forward_list.txt
ls reverse_reads/*.fq.gz > reverse_list.txt
paste forward_list.txt reverse_list.txt > file_pairs.txt
cat file_pairs.txt
```
2. Then, I will use the '--file_list' option to generate a list of input file names and check if all the file names have an R1 and R2 (pairs).  
```
fastp -i forward_reads/*.fq.gz -I reverse_reads/*.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --html report.html --json report.json
```
3. In the subdirectory containing the fastq.gz files, use the `fastp` command. Fastp is a tool that can be used for quality control and trimming of high-throughput sequencing data in fastq format, including fastq.gz compressed files.
```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --html report.html --json report.json
```
This command specifies the input files for both the forward (`-i`) and reverse (`-I`) reads, as well as the output files for the processed reads (`-o` and `-O`). The -`-html` and `--json` options specify the output files for the HTML and JSON reports, respectively.
