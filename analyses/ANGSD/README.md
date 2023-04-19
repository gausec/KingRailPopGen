## SNP filtering & variant calling
> I will be using [ANGSD (Analysis of next generation Sequencing Data)](http://www.popgen.dk/angsd/index.php/ANGSD) to do a comparative analysis of four geographically distinct populations of king rails. This will involve variant calling, genetic differentiation (FST) analysis, and estimating nucleotide diversity (π).
---
### Steps
1.  Create a VCF file containing variant calls from a BAM file using [bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html)
```
bcftools mpileup -O b -o $file_output_raw.vcf -f Rallus_crepitans_1.0.fasta -Ou --threads 20 $file.sorted.bam 
| bcftools call -m -Ov -f GQ -o $file_output_raw.vcf
```
- -f is the reference genome
- -Ou is to work with uncompressed BCF output (internal binary representation)
- -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele
- 
> *Note: I am working on this step at present. This page will be updated as I make progress.*
