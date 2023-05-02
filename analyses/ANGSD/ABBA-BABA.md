## Estimating admixture using the ABBABABA test (D-statistic) in [ANGSD](http://www.popgen.dk/angsd/index.php/Abbababa#Jackknife)

--- 


1. Make a list of the bam files (one individual from each population)
```
ls BAM/*.bam > bam_list.txt
```
---
  
  
2. ABBABABA test
```
angsd -bam bam_list.txt -doAbbababa 1 -anc CLRA.fasta -max_depth 1000 -min_maf 0.05 -block_size 5000000 -remove_bads 1 -out output
```
- `-doAbbababa`: Enables the abbababa test
- `-anc`: Specifies the reference genome file
---


3. JackKnife
```
Rscript jackKnife.R file=output.abbababa indNames=bam_list.txt outfile=out
```

- nABBA is the count of ABBA patterns, and nBABA is the count of BABA patterns.
- Dstat is a test statistic calculated (nABBA - nBABA) / (nABBA + nBABA).
  - A negative Dstat value indicates that H1 is more closely related to H3 than H2 is, and a positive value indicates that H2 is more closely related to H3 than H1 is.
- JackEst is a bias-corrected estimate of the abbababa statistic (this should be similar to the Dstat).
- SE is the estimated standard error.
- Z is the Z-score used to determine the significance of the D-statistic test. 
  - Values above 3 or below -3 are used as critical values.
