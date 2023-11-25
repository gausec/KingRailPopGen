### I am using [bcftools to calculate pairwise FST](https://speciationgenomics.github.io/variant_calling/)
---
&nbsp;

### 1. Variant Calling
&nbsp;
- 1.1 Index the Reference Genome:

```
samtools faidx OrderedCLRA.fasta
```
- 1.2. Calling Variants:

```
bcftools mpileup -a AD,DP,SP -Ou -f OrderedCLRA.fasta *.bam | bcftools call -f GQ,GP -mO z -o output.vcf.gz
```
&nbsp;
&nbsp;
### 2. Exploring VCF Files
&nbsp;
-2.1. Move VCF to a Directory:
```
mkdir vcf
mv your_output.vcf.gz ./vcf
cd vcf
```
- 2.2. Index the VCF:
```
bcftools index your_output.vcf.gz
```
&nbsp;
&nbsp;
### 3. Computing Per-Site FST
&nbsp;
- 3.1. Create a Directory for Genome Scan:
```
mkdir genome_scan
```
- 3.2. Set Environmental Variable:
```
VCF=~/vcf/your_output.vcf.gz
```
- 3.3. Extract Sample Names for Populations:
```
cd ~/genome_scan
bcftools query -l $VCF | grep "AR" > AR_pop
bcftools query -l $VCF | grep "OH" > OH_pop
bcftools query -l $VCF | grep "NC" > NC_pop
bcftools query -l $VCF | grep "FL" > FL_pop
```
- 3.4. Run vcftools for FST Calculation:
```
vcftools --gzvcf ${VCF} \
--weir-fst-pop AR_pop \
--weir-fst-pop OH_pop \
--weir-fst-pop NC_pop \
--weir-fst-pop FL_pop \
--out ./populations_fst
```
