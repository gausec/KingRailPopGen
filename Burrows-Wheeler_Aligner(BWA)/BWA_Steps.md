# Aligning 
> I am using the Burrows-Wheeler Aligner (BWA), a software tool commonly used for short read alignment to a reference genome. I will be using a [high quality clapper rail genome](https://figshare.com/articles/dataset/A_high_quality_de_novo_genome_assembly_for_Clapper_Rail_Rallus_crepitans_/21983261) (*Rallus crepitans*) to align 58 king rail (*Rallus elegans*) whole genomes. These species are sister taxa.
---
### Steps:
1. Create index files from the clapper rail genome
```
bwa index Rallus_crepitans.fasta
```
