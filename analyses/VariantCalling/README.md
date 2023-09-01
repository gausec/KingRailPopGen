## SNP filtering & variant calling
---
### Steps
1.  Create a BCF file containing variant calls from a BAM file using [bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html). 
```
bcftools mpileup -Ou -f CLRAindex/Rallus_crepitans_1.0.fasta --threads 20 11101.sorted.bam | bcftools call -mv -Ob -o "11101.bcf"
```
