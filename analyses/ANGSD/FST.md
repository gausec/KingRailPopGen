
### Aim: 
#### I will be using [ANGSD to estimate F<sub>ST</sub>](http://www.popgen.dk/angsd/index.php/Fst) between populations. [F<sub>ST</sub>](https://www.nature.com/articles/nrg2611) is a measurment of genetic variaton *within* populations reative to genetic variation *between* populations.  F<sub>ST</sub> is based on heterozygosity (a simple measure of genetic variation, 1 - the sum of the allele frequencies squared) and is reported on a scale of 0-1 (0 meaning there is 0% variation between populations and 100% variation within populations.  F<sub>ST</sub>=1  would mean there is 100% variation between and 0% among).
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
angsd -b AR_list.txt -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 1 -doSaf 1 -nThreads 10 -out FST/SFS.AR
```

&nbsp; *Repeat for the other 3 populations.*

&nbsp;

---

#### 2. Use [RealSFS](http://www.popgen.dk/angsd/index.php/RealSFS) to calculate 2d sfs for each pair of populations
```
realSFS SFS.AR.saf.idx SFS.NC.saf.idx > AR_NC.sfs; realSFS SFS.AR.saf.idx SFS.OH.saf.idx > AR_OH.sfs; realSFS SFS.AR.saf.idx SFS.FL.saf.idx > AR_FL.sfs; realSFS SFS.NC.saf.idx SFS.OH.saf.idx > NC_OH.sfs; realSFS SFS.NC.saf.idx SFS.FL.saf.idx > NC_FL.sfs; realSFS SFS.OH.saf.idx SFS.FL.saf.idx > OH_FL.sfs

```

---
#### 3. Use the above calculated 2dsfs as priors jointly with all safs from step1 to calculate fst binary files
&nbsp; 3.1
```
```
- *I am using the* `-fold 1` *option for a single SAF file that includes all populations.*

---
#### 4. Use realSFS to extract the the fst values from the fst binary files
&nbsp; 4.1


