# Aligning to reference genome
> I am using the Burrows-Wheeler Aligner (BWA), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://figshare.com/articles/dataset/A_high_quality_de_novo_genome_assembly_for_Clapper_Rail_Rallus_crepitans_/21983261) (*Rallus crepitans*) to align 58 king rail (*Rallus elegans*) whole genomes. These species are sister taxa.
---
### Steps:
1. Create index directory & index files from the clapper rail genome
```
mkdir index
bwa index index/Rallus_crepitans.fasta
```

2. Assign basename & align the trimmed reads to the reference genome using the read alignment tool
```
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 16 index/Rallus_crepitans.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > SAM/${base}.sam; done &
```
*Note: Above, the loop iterates through each R1 fastq file in the current directory using the wildcard `*R1.fastq.gz`. For each R1 file, the base filename is extracted using the `basename` command, and the `-R1.fastq.gz` is removed using the `-s` option. This base name is then used to make the names of the corresponding R2 file and the output SAM file. Finally, the BWA alignment is performed using `bwa mem` with the clapper rail reference genome & input fastq files. The output is redirected to a SAM file using the `>` operator. The ampersand (&) at the end runs the command in the background.*


3. Convert SAM to BAM using samtools (this saves a lot of space by converting to binary data)
```
for r1file in *.sam; do base=$(basename $r1file .sam); samtools view -@ 40 $r1file > ${base}.bam; done
```
