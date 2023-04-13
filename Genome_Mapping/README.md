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
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 16 index/Rallus_crepitans_1.0.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done &
```
- *Note: Above, the loop iterates through each R1 fastq file in the current directory using the wildcard `*R1.fastq.gz`. For each R1 file, the base filename is extracted using the `basename` command, and the `-R1.fastq.gz` is removed using the `-s` option. This base name is then used to make the names of the corresponding R2 file and the output SAM file.*
- *Finally, the BWA alignment is performed using `bwa mem` with the clapper rail reference genome & input fastq files. To speed up the process, `-t` indicates I want to use 16 threads (each thread is assigned a portion of the total computing resources (like CPU time, memory) and can execute independently of each other). The output is redirected to a SAM file using the `>` operator. The ampersand (&) at the end runs the command in the background.*
- *Megan: don't worry about the `&` since you have 4  files. I just did this so I can disown the process and leave it to run.*

---
#### 3. Convert SAM to BAM using [samtools](https://github.com/samtools/samtools) (this saves a lot of space by converting to binary data)
```
for r1file in *.sam; do base=$(basename $r1file .sam); samtools view -@ 40 $r1file > ${base}.bam; done
```
---
#### 4. Arrange the aligned reads in the BAM files based on their coordinates in the reference genome (This is called sorting)
```
for file in 11101.bam 11103.bam 11104.bam ... RK02.bam; do echo "Sorting ${file}..."; samtools sort -@ 40 -T temp_dir ${file} -o sorted/${file%.*}.sorted.bam; done
```
*Note: The ellipsis indicates that there are similar file types. I am using `echo` so bash will tell me which file samtools is currently sorting. `-@` specifies that I want to use 40 threads to speed up the process. `${file%.*}` is a parameter expansion that removes the shortest matching pattern from the variable `$file`. In this case, the pattern is `.*`. This matches any sequence of characters that ends with a dot and removes it from the end of `$file`.* 
 
#### 5. Determine average depth
```
samtools depth 11101.sorted.bam  11111.sorted.bam  11120.sorted.bam  43633.sorted.bam  43670.sorted.bam  49805.sorted.bam  6015.sorted.bam  6026.sorted.bam 11103.sorted.bam  11113.sorted.bam  43215.sorted.bam    43639.sorted.bam  43671.sorted.bam  49806.sorted.bam  6016.sorted.bam  6027.sorted.bam 11104.sorted.bam  11114.sorted.bam  43613.sorted.bam    43646.sorted.bam  45049.sorted.bam  50753.sorted.bam  6017.sorted.bam  RK01.sorted.bam 11105.sorted.bam 11115.sorted.bam  43617.sorted.bam 43647.sorted.bam 45309.sorted.bam  50930.sorted.bam 6018.sorted.bam  RK02.sorted.bam 11108.sorted.bam 11116.sorted.bam  43621.sorted.bam 43649.sorted.bam  48819.sorted.bam  51747.sorted.bam 6019.sorted.bam 11109.sorted.bam 11117.sorted.bam  43623.sorted.bam 43654.sorted.bam 48828.sorted.bam  53078.sorted.bam 6020.sorted.bam 11110.sorted.bam 11119.sorted.bam 43657.sorted.bam 49804.sorted.bam 6014.sorted.bam 6025.sorted.bam | awk '{sum+=$3} END {print "Average depth:", sum/NR}' > AverageDepth.txt
```
