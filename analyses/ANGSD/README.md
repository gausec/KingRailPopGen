## SNP filtering & variant calling
> I will be using [ANGSD (Analysis of next generation Sequencing Data)](http://www.popgen.dk/angsd/index.php/ANGSD) to do a comparative analysis of four geographically distinct populations of king rails. This will involve variant calling, genetic differentiation (FST) analysis, and estimating nucleotide diversity (π).
---
### Steps
1. Using
```
bcftools mpileup -f Rallus_crepitans_1.0.fasta -Ou -q 20 -Q 20 -C 50 -I -d 1000000 -a AD $file.sorted.bam | bcftools call -mv -Ov -f GQ -o $file_output_raw.vcf
```
