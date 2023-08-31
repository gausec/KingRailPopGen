## SNP filtering & variant calling
---
### Steps
1.  Create a VCF file containing variant calls from a BAM file using [bcftools](https://samtools.github.io/bcftools/howtos/variant-calling.html). 
```
./bcftools.sh
```
- *This script can be found in the bin subdirectory*
2. bcftools can be used to query a vcf file. Here I am searching for the B-95207 isolate:
```
bcftools query -i 'ID="CM053018.1"' -f '%CHROM:%POS\n' 11101_output_raw.vcf > output.txt
```
