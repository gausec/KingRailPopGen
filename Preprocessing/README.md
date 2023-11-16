$${\color{green}\Huge {Preprocessing \space for \space FastQ \space files }}$$   

> This subdirectory contains steps for the preprocessing of paired-end whole genome data of 4 species (king rail, clapper rail, Virginia rail, and eastern bluebird). All steps below occured in the command line once the original WGS data was saved to an external hard drive. The forward and reverse fastq files were in two separate tar files.
---
### Step 1: File transfer, extraction, and organization

1.1 Create a subdirectory within the remote server
```
mkdir WholeGenomes 
```
1.2 Transfer tar files to remote server using the `scp` command

1.3 Extract folders containing zipped fastq files 
```
tar -xvf sequences.tar
```
1.4 Navigate to these R1 & R2 folders using the `cd` command, and check that the fastq.gz files were successfully extracted

1.5 Navigate back to previous directory using `cd ..` command and remove the tar files to save space
```
rm sequences.tar
```
1.6 Create a subdirectory for the cleaned reads and another for the json and html reports in both R1 & R2 folders
```
mkdir cleaned
mkdir reports
```
*Note: fastp was installed this project. [fastp]( https://github.com/OpenGene/fastp) is a tool that can be used for quality control and trimming of high-throughput sequencing data in fastq format, including fastq.gz compressed files.*

---
### Step 2: Quality control

2.1 Rename files to simplify downstream tasks. Use the `cut` command to extract the first 5 characters of the file name, and store it as the variable `name`.
```
for FileName in *_R1_001.fastq.gz; do name=$(echo "$FileName" | cut -c1-5); mv "$FileName" "${name}_R1.fastq.gz"; done
```
```
for FileName in *_R2_001.fastq.gz; do name=$(echo "$FileName" | cut -c1-5); mv "$FileName" "${name}_R2.fastq.gz"; done  
```
2.2 Use a for loop to run the `fastp` command (without changing default quality parameters) for all fastq files. This command specifies the input files for both the forward (`-i`) and reverse (`-I`) reads, as well as the output files for the processed reads (`-o` and `-O`).  

```
for FileName in *_R1.fastq.gz; do name=$(echo $FileName | cut -c1-5); fastp -i $name*_R1.fastq.gz -I $name*_R2.fastq.gz -o cleaned/$name*-clean-R1.fastq.gz -O cleaned/$name-clean-R2.fastq.gz -w 10 --html reports/$name.html --json reports/$name.json; done
```
*Note: The above code contains a for loop that uses the variable `FileName` to iterate over all files ending in `_R1_001.fastq.gz` in the current directory. For each file, it extracts the first part of the file name (before the first underscore) using the `ls` and `cut` commands. This is saved to the variable name. Then, it runs the fastp command using the input file names with the `$name` variable inserted in appropriate positions, and puts the cleaned files in the `/cleaned` directory. The `--html` and `--json` options are used to generate an HTML and a JSON report file for each input file.* 


