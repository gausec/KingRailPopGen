# Aligning to reference genome
> I am using the Burrows-Wheeler Aligner (BWA), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://figshare.com/articles/dataset/A_high_quality_de_novo_genome_assembly_for_Clapper_Rail_Rallus_crepitans_/21983261) (*Rallus crepitans*) to align 58 king rail (*Rallus elegans*) whole genomes. These species are sister taxa.
---
### Steps:
1. Create index files from the clapper rail genome
```
bwa index Rallus_crepitans.fasta
```
2. Assign basename & align the trimmed reads to the reference genome using the read alignment tool
```
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 16 index/Rallus_crepitans.fasta ${base}-R1.fastq.gz ${base}-R2.fastq.gz > ${base}.sam; done
```
