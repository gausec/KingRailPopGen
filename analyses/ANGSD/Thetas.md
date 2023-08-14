### I will be estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests).
---

> Note: This method requires a SAF file for a given population. I created an SAF index file for each population when I estimated [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md). **All steps below were repeated for each population.**
&nbsp;
&nbsp;
---
&nbsp;
#### 1.Estimate the unfolded site frequency spectrum
```
../../../../angsd/misc/realSFS SFS.NC.saf.idx -P 24 > out.NC.sfs
```
&nbsp;

#### 2. Calculate per-site thetas
```
../../../../angsd/misc/realSFS saf2theta SFS.NC.saf.idx -P 20 -sfs out.NC.sfs -outname NC.out
```

&nbsp;
#### 3. Calculate neutrality tests statistics

```
../../../../angsd/misc/thetaStat do_stat NC.out.thetas.idx
```
&nbsp;
#### 4. Use the [thetaStat program](http://www.popgen.dk/angsd/index.php/ThetaStat) to estimate Tajima’s D.
```

```
---
&nbsp;
### Next, I want to calculate nucleotide diversity ($\pi$) following the same method as [this study](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-018-1209-y).

&nbsp;
#### 5. Extract the theta P (tP) column from the *.thetas.gz.pestPG output file and divide by the number of sites (nSites) used for the population.
```
awk '{print $5 / $14}' NC.out.thetas.idx.pestPG > tP_divided_by_numSites.txt

```

