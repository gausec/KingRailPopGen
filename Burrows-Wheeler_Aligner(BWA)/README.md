# Aligning to reference genome
> I am using the Burrows-Wheeler Aligner (BWA), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://figshare.com/articles/dataset/A_high_quality_de_novo_genome_assembly_for_Clapper_Rail_Rallus_crepitans_/21983261) (*Rallus crepitans*) to align 58 king rail (*Rallus elegans*) whole genomes. These species are sister taxa.


### Steps:

---
1. Create index directory & index files from the clapper rail genome
```
mkdir index
bwa index index/Rallus_crepitans.fasta
```
---
2. Assign basename & align the trimmed reads to the reference genome using the read alignment tool
```
mkdir SAM
```
```
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 16 index/Rallus_crepitans.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done &
```
*Note: Above, the loop iterates through each R1 fastq file in the current directory using the wildcard `*R1.fastq.gz`. For each R1 file, the base filename is extracted using the `basename` command, and the `-R1.fastq.gz` is removed using the `-s` option. This base name is then used to make the names of the corresponding R2 file and the output SAM file. Finally, the BWA alignment is performed using `bwa mem` with the clapper rail reference genome & input fastq files. The output is redirected to a SAM file using the `>` operator. The ampersand (&) at the end runs the command in the background-- Megan: don't worry about the `&` since you have 4  files. I just do this so I can disown the process and leave it to run.*

---
3. Convert SAM to BAM using samtools (this saves a lot of space by converting to binary data)
```
for r1file in *.sam; do base=$(basename $r1file .sam); samtools view -@ 40 $r1file > ${base}.bam; done
```
---
4. Arrange the aligned reads in the BAM files based on their coordinates in the reference genome (This is called sorting)
```
for file in 11101.bam 11103.bam 11104.bam ... RK02.bam; do echo "Sorting ${file}..."; samtools sort -@ 40 -T temp_dir ${file} -o sorted/${file%.*}.sorted.bam; done
```
*Note: The ellipsis indicates that there are similar file types. I am using `echo` so bash will tell me which file samtools is currently sorting. `-@` specifies that I want to use 40 threads to speed up the process. `${file%.*}` is a parameter expansion that removes the shortest matching pattern from the variable `$file`. In this case, the pattern is `.*`. This matches any sequence of characters that ends with a dot and removes it from the end of `$file`.*
