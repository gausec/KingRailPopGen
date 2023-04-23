# Mapping Reads to a Reference Genome
> I am using the [Burrows-Wheeler Aligner (BWA)](https://github.com/lh3/bwa), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://figshare.com/articles/dataset/A_high_quality_de_novo_genome_assembly_for_Clapper_Rail_Rallus_crepitans_/21983261) (*Rallus crepitans*) to align 57 king rail (*Rallus elegans*) whole genomes. These species are sister taxa. 


## Steps:

#### 1. Create index directory & index files from the clapper rail genome
```
mkdir index
bwa index index/Rallus_crepitans_1.0.fasta
```
---
#### 2. Assign basename & align the trimmed reads to the reference genome using the read alignment tool
```
mkdir SAM
```
```
screen # running the code on a background screen
```
```
for file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${file}); bwa mem -t 16 index/Rallus_crepitans_1.0.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done
```
- Next press ctrl + `a` and then `d` to disown the process and then leave close the terminal to let the command finish running.
- *Note: Above, the loop iterates through each R1 fastq file in the current directory using the wildcard `*R1.fastq.gz`. For each R1 file, the base filename is extracted using the `basename` command, and the `-R1.fastq.gz` is removed using the `-s` option. This base name is then used to make the names of the corresponding R2 file and the output SAM file.*
- *Finally, the BWA alignment is performed using `bwa mem` with the clapper rail reference genome & input fastq files. To speed up the process, `-t` indicates I want to use 16 threads (each thread is assigned a portion of the total computing resources (like CPU time, memory) and can execute independently of each other). The output is redirected to a SAM file using the `>` operator. The ampersand (&) at the end runs the command in the background.*

---
#### 3. Convert SAM to BAM using [samtools](https://github.com/samtools/samtools) (this saves a lot of space by converting to binary data)
```
for file in *.sam; do base=$(basename $file.sam); samtools view -@ 40 $file > ${base}.bam; done
```

---
#### 4. Arrange the aligned reads in the BAM files based on their coordinates in the reference genome (This is called sorting)
```
for file in *.bam; do echo "Sorting ${file}"; samtools sort -@ 40 -T temp_dir ${file} -o sorted/${file%.*}.sorted.bam; done
```
*Note: I am using `echo` so bash will tell me which file samtools is currently sorting. `-@` specifies that I want to use 40 threads to speed up the process. `${file%.*}` is a parameter expansion that removes the shortest matching pattern from the variable `$file`. In this case, the pattern is `.*`. This matches any sequence of characters that ends with a dot and removes it from the end of `$file`.* 
 
#### 5. Determine average depth
```
for file in *.sorted.bam; do echo "$file" >> AvgDepth.txt; samtools depth $file | awk '{sum+=$3} END {print "Average depth:", sum/NR}' >> AvgDepth.txt;  echo "" >> AvgDepth.txt; done
```
*This loop prints the filename to a text file named AvgDepth.txt, runs samtools depth on the file & pipes the output to 'awk' to calculate the average depth and outputs it to the same text file. It then prints an empty line to separate the text. The loop ends when all the .bam files have been processed.* 

#### 6. Mapping summary stats
```
for file in *.sorted.bam; do echo "$file" >> KIRAflagstats.txt; samtools flagstat -@ 30 "$file" >> KIRAflagstats.txt; echo "" >> KIRAflagstats.txt; done
```
*This loop prints the filename to a text file named KIRAflagstats.txt, runs samtools flagstat on the file & appends the output to the same text file, and then prints an empty line to separate the text. The loop ends when all the .bam files have been processed.*
