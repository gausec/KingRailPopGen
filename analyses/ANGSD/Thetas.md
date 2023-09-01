## Estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
---
Note: This method requires a SFS file for a given population. Tajima’s D has to be folded, so the SFS index file for each population that I created when estimating [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md) can't be used here. **All steps below were repeated for each population.**
&nbsp;

--- 
&nbsp;

#### 1. Estimate the folded site frequency spectrum.
```
../../../../angsd/misc/realSFS SFS.NC.saf.idx -P 24 -fold 1 -anc ../../CLRAindex/Rallus_crepitans_1.0.fasta > out.NC.sfs
```
&nbsp;

#### 2. Calculate per-site thetas.
```
../../../../angsd/misc/realSFS saf2theta SFS.NC.saf.idx -P 20 -sfs out.NC.sfs -outname NC.out
```

&nbsp;
#### 3. Calculate neutrality test statistics using the [thetaStat program](http://www.popgen.dk/angsd/index.php/ThetaStat).
```
../../../../angsd/misc/thetaStat do_stat NC.out.thetas.idx
```
&nbsp;
#### 4. Sliding window analysis to calculate statistics in genomic windows.
```
 ../../../../angsd/misc/thetaStat do_stat NC.out.thetas.idx -win 5000 -step 1000 -outnames NC
```
- The *.thetas.idx.pestPG file contains statistics for different genomic regions.
- Columns contain information about regions, reference name, center of the window, estimators of theta (Watterson, pairwise, etc.), and neutrality test statistics (Tajima's D, Fu&Li F's, etc.).
- The last column indicates the effective number of sites with data in the window.
---
&nbsp;
### Next, I want to calculate nucleotide diversity ($\pi$) following the same method as [this study](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-018-1209-y), and then calculate Watterson's theta.

&nbsp;

#### 5. Extract the theta P (tP) column from the *.thetas.gz.pestPG output file and divide by the number of sites (nSites) used for the population.
```
awk '{print $5 / $14}' NC.out.thetas.idx.pestPG >> pi.txt
```

#### 6. Remove any empty rows.
```
awk '$1 != "-nan"' pi.txt > cleaned_pi.txt
```

#### 7. Calculate pi.
```
awk '{ sum += $1 } END { print sum / NR }' cleaned_pi.txt
```
#### 8. Calculate Watterson's theta from the *.thetas.gz.pestPG output file in the same way.
8.1 Extract Watterson's theta (tW).
```
awk '{print $5 / $14}' NC.out.thetas.idx.pestPG >> tW.txt
```
8.2 Remove any empty rows.
```
awk '$1 != "-nan"' tW.txt > cleaned_tW.txt
```
8.3 Calculate Waterson's theta by taking the average.
```
awk '{ sum += $1 } END { print sum / NR }' cleaned_tW.txt
```
