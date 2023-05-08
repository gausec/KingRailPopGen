# Fst

---

### I will be using ANGSD to calculated [Fst](http://www.popgen.dk/angsd/index.php/Fst)
---
### 1. SFS estimation for each population
1.1 Create a text file list of bam files for each population. These can be found in the data subdirectory. 
- NC_list.txt
- OH_list.txt
- AR_list.txt
- FL_list.txt

1.2 Generage site allele frequency likelihoods in SAF format for each population.
```
angsd -b AR_list.txt -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 1 -doSaf 1 -out FST/SFS.AR
```
*Repeat for the other 3 populations.*
