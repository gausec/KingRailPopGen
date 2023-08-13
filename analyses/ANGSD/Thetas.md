### I will be estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests).
---

> Note: This method requires a SAF file for a given population. I created an SAF index file for each population when I estimated [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md). **All steps below were repeated for each population.**
&nbsp;
&nbsp;
---
&nbsp;
### 1.Estimate the unfolded site frequency spectrum
```
../../../../angsd/misc/realSFS SFS.NC.saf.idx -P 24 > out.NC.sfs
```
&nbsp;

### 2. Calculate per-site thetas
```
../../../../angsd/misc/realSFS saf2theta SFS.NC.saf.idx -P 20 -sfs out.NC.sfs -outname NC.out
```

&nbsp;
### 3. Estimate Tajimas D and other statistics

```
./misc/thetaStat do_stat out.thetas.idx
 ```

