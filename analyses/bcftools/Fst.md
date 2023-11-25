## 1. Variant Calling
- Index the Reference Genome:

```
samtools faidx ImprovedCLRA.fasta
```
- 1.2. Calling Variants:

```
REF=~/reference/ImprovedCLRA.fasta
bcftools mpileup -a AD,DP,SP -Ou -f $REF ${ARlist[@]} ${OHlist[@]} ${NClist[@]} ${FLlist[@]} | bcftools call -f GQ,GP -mO z -o ./your_output.vcf.gz
```

## 2. Exploring VCF Files:
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
## 3. Computing Per-Site FST:
- 3.1. Create a Directory for Genome Scan:
```
mkdir genome_scan
```
-3.2. Set Environmental Variable:
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
