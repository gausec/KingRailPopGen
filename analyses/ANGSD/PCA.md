## Principal component analysis
#### I am using principal component analysis through [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd) to better understand population structure.
---
&nbsp;
### Steps: 
1. Prepare input data by obtaining genotype likelihoods in Beagle format from BAM files
```
../../angsd/angsd -GL 2 -out ../../PCA/AR_genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam AR_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/NC_genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam NC_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/FL_genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam FL_list.txt;
../../angsd/angsd -GL 2 -out ../../PCA/OH_genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam OH_list.txt
```
2. Run PCAngsd
```
python pcangsd.py -beagle ../../PCA/input.beagle.gz -out output -threads 10
```
