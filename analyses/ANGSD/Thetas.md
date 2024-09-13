## Estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
---
Note: This method requires a SFS file for a given population. Tajima’s D has to be folded, so the SFS index file for each population that I created when estimating [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md) can't be used here. **All steps below were repeated for each population.**
&nbsp;

--- 
&nbsp;

#### 1. Estimate the folded site frequency spectrum.
&nbsp; 1.1 Generate site allele frequency likelihoods in SAF format for each population.
```
angsd/angsd -b OHlist -anc Ordered.CLRA.fasta -ref Ordered.CLRA.fasta -minMapQ 30 -minQ 20 -GL 1 -minInd 8 -doSaf 1 -baq 1 -nthreads 8 -out SAF.OH 
```
&nbsp; *Repeat for the other 3 populations.*

- `-GL 1`: SAMtools model: GATK assumes sequencing errors are independent, while Samtools believes the second error (unrelated sequencing error that occurs independently of the first one and affects the same genomic position) comes at a higher chance, especially at high depths.

&nbsp;

&nbsp; 1.2 Generate folded SFS for each population.
```
angsd/misc/realSFS SAF.OH.saf.idx -P 24 -fold 1 -anc Ordered.CLRA.fasta -ref Ordered.CLRA.fasta > Ohio.Folded.SFS 
```
&nbsp;

#### 2. Calculate per-site thetas.
```
angsd/misc/realSFS saf2theta SAF.OH.saf.idx -outname Ohio.Folded -sfs Ohio.Folded.SFS -fold 1
```

&nbsp;
#### 3. Calculate neutrality test statistics using the [thetaStat program](http://www.popgen.dk/angsd/index.php/ThetaStat).
```
angsd/misc/thetaStat do_stat Ohio.Folded.thetas.idx 
```
&nbsp;
#### 4. Sliding window analysis to calculate statistics in genomic windows. I skipped this step.
```
angsd/misc/thetaStat do_stat NC.out.thetas.idx -win 5000 -step 1000 -outnames NC
```
- The *.thetas.idx.pestPG file contains statistics for different genomic regions.
- Columns contain information about regions, reference name, center of the window, estimators of theta (Watterson, pairwise, etc.), and neutrality test statistics (Tajima's D, Fu&Li F's, etc.).
- The last column indicates the effective number of sites with data in the window.
&nbsp;
&nbsp;

---

$${\color{orange}Nucleotide \space diversity \space (π) \space and \space Watterson's \space theta}$$

---
### Next, I want to calculate nucleotide diversity (π) following the same method as [this study](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-018-1209-y), and then calculate Watterson's theta. *There are now scripts for these in the bin subdirectory.*

&nbsp;

#### 5. Extract the theta P (tP) column from the *.thetas.gz.pestPG output file and divide by the number of sites (nSites) used for the population. 
```
awk 'NR > 1 {print $5 / $14}' NC.pestPG >> pi.txt
```

#### 6. Calculate average pi.
```
awk '{ sum += $1 } END { print sum / NR }' pi.txt
```
#### 7. Calculate Watterson's theta from the *.thetas.gz.pestPG output file in the same way.
7.1 Extract Watterson's theta (tW).
```
awk 'NR > 1 {print $4 / $14}' NC.pestPG >> tW.txt
```
7.2 Calculate Waterson's theta by taking the average.
```
awk '{ sum += $1 } END { print sum / NR }' tW.txt
```
 
