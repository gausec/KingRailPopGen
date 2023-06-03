
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
../../angsd/angsd -b NC_list.txt -anc CLRAindex/Rallus_crepitans_1.0.fasta -gl 2 -doSaf 1 -doMaf 1 -SNP_pval 0.05 -minind 5 -nthreads 8 -doMajorMinor 1 -out FST/SFS.NC

```

&nbsp; *Repeat for the other 3 populations.*

&nbsp;

---

#### 2. Use [RealSFS](http://www.popgen.dk/angsd/index.php/RealSFS) to calculate 2d sfs for each pair of populations
```
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx NC/SFS.NC.saf.idx -P 30 > AR_NC.ml;
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx OH/SFS.OH.saf.idx -P 30 > AR_OH.ml; 
../../../angsd/misc/realSFS AR/SFS.gl2.AR.saf.idx FL/SFS.FL.saf.idx -P 30 > AR_FL.ml; 
../../../angsd/misc/realSFS NC/SFS.NC.saf.idx OH/SFS.OH.saf.idx -P 30 > NC_OH.ml; 
../../../angsd/misc/realSFS NC/SFS.NC.saf.idx FL/SFS.FL.saf.idx -P 30 > NC_FL.ml; 
../../../angsd/misc/realSFS OH/SFS.OH.saf.idx FL/SFS.FL.saf.idx -P 30 > OH_FL.ml

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


