# Mapping Reads to a Reference Genome
> I am using the [Burrows-Wheeler Aligner (BWA)](https://github.com/lh3/bwa), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://www.ncbi.nlm.nih.gov/assembly/GCA_028554615.1/) (*Rallus crepitans*) to align 57 king rail (*Rallus elegans*) whole genomes. These species are sister taxa. 


## Steps:

#### 1. Create index directory & index files from the clapper rail genome
```
mkdir index
bwa index index/Rallus_crepitans_1.0.fasta
```
---

&nbsp;

#### 2. Assign basename & align the trimmed reads to the reference genome using the read alignment tool
```
mkdir SAM
```
2.1 Run the code in the background. The base name is extracted and stored as a local variable, `base`
```
screen
```
```
for file in *R1.fastq.gz; do base=$(basename "$file" | cut -c1-5); echo "$base"; bwa mem -t 16 index/Rallus_crepitans_1.0.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done
```
 -  *Next use `ctrl a` and then `d` to disown the process. The terminal can now be closed, and the command will finish running.*
 -  *BWA alignment is performed using `bwa mem` with the clapper rail reference genome & fastq files. To speed up the process, `-t` indicates I want to use 16 threads (each thread is assigned a portion of the total computing resources, such as CPU time, and can execute independently of each other). The output is redirected to a SAM file using the `>` operator.*

---

&nbsp;

#### 3. Convert SAM to BAM using [samtools](https://github.com/samtools/samtools) (this saves a lot of space by converting to binary data)
```
for file in *.sam; do base=$(basename "$file" .sam); samtools view -bS -@ 8 $file > BAM/${base}.bam; done
```
- *Note: to save time, pipe the mapped output directly to samtools to sort and save output to a bam file.*
- *Below is an example using the improved reference genome post-RagTag*:
```
bwa mem -t 30 OrderedCLRA.fasta 11101-clean-R1.fastq.gz 11101-clean-R2.fastq.gz | samtools sort -@ 20 -o 11101.Improved.bam -
```

---

&nbsp;

#### 4. Arrange the aligned reads in the BAM files based on their coordinates in the reference genome (This is called sorting)
```
for file in *.bam; do echo "Sorting ${file}"; samtools sort -@ 8 -o SortedBAM/${file}.sorted.bam; done
```
*Note: I am using `echo` so bash will tell me which file samtools is currently sorting. `-@` specifies that I want to use threads. `${file%.*}` is a parameter expansion that removes the shortest matching pattern from the variable `$file`. In this case, the pattern is `.*`. This matches any sequence of characters that ends with a dot and removes it from the end of `$file`.* 

 &nbsp;
 
#### 5. Determine average depth
```
for file in *.sorted.bam; do echo "$file" >> AvgDepth.txt; samtools depth $file | awk '{sum+=$3} END {print "Average depth:", sum/NR}' >> AvgDepth.txt;  echo "" >> AvgDepth.txt; done
```
*This loop prints the filename to a text file named AvgDepth.txt, runs samtools depth on the file & pipes the output to 'awk' to calculate the average depth and outputs it to the same text file. It then prints an empty line to separate the text. The loop ends when all the .bam files have been processed.* 

&nbsp;

#### 6. Mapping summary stats
```
for file in *.sorted.bam; do echo "$file" >> KIRAflagstats.txt; samtools flagstat -@ 30 "$file" >> KIRAflagstats.txt; echo "" >> KIRAflagstats.txt; done
```
*This loop prints the filename to a text file named KIRAflagstats.txt, runs samtools flagstat on the file & appends the output to the same text file, and then prints an empty line to separate the text. The loop ends when all the .bam files have been processed.*

&nbsp;

#### 3. Sanity check using [samtools quickcheck](https://www.htslib.org/doc/samtools-quickcheck.html)
```
samtools quickcheck -v *.bam > bad_bams.fofn   && echo 'all files passed the check' || echo 'some files failed the check, see bad_bams.fofn'
```
Notes about sanity check script:
- The `-v` flag enables verbose output, which provides more details about any issues.
- The `.fofn` extension stands for "File of File Names." It's a convention used to indicate that the file contains a list of file names, often one per line. Here, `bad_bams.fofn` will contain a list of any files that failed the samtools quickcheck.
- The double pipe `||` in shell scripting runs the second command only if the first command fails. It’s a way to handle errors or run an alternative action if the first command doesn’t work.
&nbsp;
