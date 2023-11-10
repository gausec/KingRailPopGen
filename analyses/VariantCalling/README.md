## SNP calling using bcftools
---
### Steps
1.  Create a BCF file containing variant calls from a BAM file using [bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html). 
```
bcftools mpileup -Ou -f CLRAindex/Rallus_crepitans_1.0.fasta --threads 20 11101.sorted.bam | bcftools call -mv -Ob -o "11101.bcf"
```

---

### SNP calling using ANGSD
```
angsd/angsd -bam BamFileList -ref CLRAindex/Rallus_crepitans_1.0.fasta -anc CLRAindex/Rallus_crepitans_1.0.fasta -GL 2 -doMajorMinor 1 -doMaf 2 -doSaf 1 -minMaf 0.05 -minind 5 -SNP_pval 1e-6 -minMapQ 30 -minQ 20 -out input -nThreads 8
```

- `-bam BamFileList` = analyze data from bam files 
- `-GL 2` = calculate the genotype likelihood using the GATK method (`-GL 1` would be SAMtools model)
- `-doMajorMinor 1` = infer the major and minor alleles
- `-doMAF 2` = estimate the allele frequencies assuming known minor 
- `-SNP_pval 1e-6` = only keep those sites that have a p-value less than 1e-6 for being variable
