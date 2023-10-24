## NGSadmix in ANGSD

1. Create the the beagle genotype likelihood input file using ANGSD: calculate genotype likelihoods for polymorphic sites using ANGSD (NGSadmix uses Genotype Likelihoods (GLs) in .beagle format as input)
```
ANGSD=~/angsd/angsd
```
```
$ANGSD -bam all.files -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -doGlf 2 -out input.gz -P 5
```
2. Run an analysis of the GLs with NGSadmix, assuming the number of ancestral populations is K:
3. ```
   $NGSADMIX -likes input.gz -K 1 -minMaf 0.05 -seed 1 -o $OUT_DIR/NGSadmix
   ```
