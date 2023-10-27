## NGSadmix in ANGSD

1. Create the the beagle genotype likelihood input file using ANGSD: calculate genotype likelihoods for polymorphic sites using ANGSD (NGSadmix uses Genotype Likelihoods (GLs) in .beagle format as input)
```
ANGSD=~/angsd/angsd
```
```
$ANGSD -bam BamFileList -GL 1 -doMajorMinor 1 -doMaf 1 -ref CLRAindex/Rallus_crepitans_1.0.fasta -SNP_pval 1e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -doGlf 2 -out input.gz -nThreads 8
```
- -GL 1 = SAMtools model
- -doGlf 3 = beagle binary
2. Run an analysis of the GLs with NGSadmix, assuming the number of ancestral populations is K:
```
 for ((k=1; k<=10; k++)); do ../../../angsd/NGSadmix -likes input.gz.glf.gz -K $k -minMaf 0.05 -seed 1 -o NGSadmix_$k; done
```
- Trying `-k` = 1-10
## Plot the results in R
