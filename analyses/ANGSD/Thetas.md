## Estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
---
Note: This method requires a SFS file for a given population. Tajima’s D has to be folded, so the SFS index file for each population that I created when estimating [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md) can't be used here. **All steps below were repeated for each population.**
&nbsp;

--- 
&nbsp;

#### 1. Estimate the folded site frequency spectrum
```
../../../../angsd/misc/realSFS SFS.NC.saf.idx -P 24 -fold 1 -anc ../../CLRAindex/Rallus_crepitans_1.0.fasta > out.NC.sfs
```
&nbsp;

#### 2. Calculate per-site thetas
```
../../../../angsd/misc/realSFS saf2theta SFS.NC.saf.idx -P 20 -sfs out.NC.sfs -outname NC.out
```

&nbsp;
#### 3. Calculate neutrality tests statistics using the [thetaStat program](http://www.popgen.dk/angsd/index.php/ThetaStat).
```
../../../../angsd/misc/thetaStat do_stat NC.out.thetas.idx
```
&nbsp;
#### 4. Sliding window analysis
```
 ../../../../angsd/misc/thetaStat do_stat FL.out.thetas.idx -win 5000 -step 1000 -outnames NC
```
---
&nbsp;
### Next, I want to calculate nucleotide diversity ($\pi$) following the same method as [this study](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-018-1209-y).

&nbsp;
#### 5. Extract the theta P (tP) column from the *.thetas.gz.pestPG output file and divide by the number of sites (nSites) used for the population
```
awk '{print $5 / $14}' NC.out.thetas.idx.pestPG > Pi.txt

```

