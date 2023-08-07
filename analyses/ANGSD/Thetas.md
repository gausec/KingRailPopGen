### I will be estimating Thetas values and other neutrality statistics using [ANGSD](http://popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests).

> Note: This method requires output from the site frequency spectrum to run. I calculated the SFS for each population when I estimated [Fst](https://github.com/gausec/KingRailPopGen/blob/main/analyses/ANGSD/FST.md).
>
> ### 1. Calculate per-site thetas
> ```
> realSFS saf2theta out.saf.idx -sfs out.sfs -outname out
> ```
>
> ### 2. Estimate Tajimas D and other statistics
> ```
> ./misc/thetaStat do_stat out.thetas.idx
> ```
