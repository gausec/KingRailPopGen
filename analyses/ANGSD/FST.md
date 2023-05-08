
### Aim: 
#### I will be using [ANGSD to estimate F<sub>ST</sub>](http://www.popgen.dk/angsd/index.php/Fst) between populations. [F<sub>ST</sub>](https://www.nature.com/articles/nrg2611) is a measurment of genetic variaton within populations reative to genetic variation between populations.  F<sub>ST</sub> is based on heterozygosity (a simple measure of genetic variation, 1 - the sum of the allele frequencies squared) and is reported on a scale of 0-1 (0 meaning there is 0% variation between populations and 100% variation within populations.  F<sub>ST</sub>=1  would mean there is 100% variation between and 0% among).
&nbsp;
---
### Steps:  
#### 1. SFS estimation for each population
&nbsp; 1.1 Create a text file list of bam files for each population. 
- NC_list.txt
- OH_list.txt
- AR_list.txt
- FL_list.txt  

&nbsp; *These can be found in the data subdirectory.*
      
&nbsp;

&nbsp; 1.2 Generate site allele frequency likelihoods in SAF format for each population.
```
angsd -b AR_list.txt -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 1 -doSaf 1 -out FST/SFS.AR
&nbsp;

```
&nbsp; *Repeat for the other 3 populations.*

---
#### 2.
&nbsp; 2.1 
