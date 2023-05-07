# Fst

---

### I will be using ANGSD to calculated Fst
---
### 1. SFS estimation for each population
1.1 Create a text file list of bam files for each population


```
angsd -b bam_list.txt -anc CLRA.fasta -gl 1 -doSaf 1 -out SFS
```
