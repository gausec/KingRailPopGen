## NGSadmix
&nbsp;
> I am using [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmixTutorial) in ANGSD to estimate individual admixture proportions from NGS data for  four populations of king rails, and then examine whether these populations separate based on their admixture patterns. I am looking for substructure among the populations.

&nbsp;

Helpful tutorial: [Population Structure using NGSadmix](https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html)

&nbsp;
&nbsp;
#### 1. Create the the beagle genotype likelihood input file using ANGSD: calculate genotype likelihoods for polymorphic sites using ANGSD (NGSadmix uses Genotype Likelihoods (GLs) in .beagle format as input)
1.1 shortcut
```
ANGSD=../../angsd
```
1.2 Beagle file
```
$ANGSD/angsd -bam BamFileList -ref CLRAindex/Rallus_crepitans_1.0.fasta -anc CLRAindex/Rallus_crepitans_1.0.fasta -GL 1 -doMajorMinor 1 -doMaf 2 -minMaf 0.05 -minind 5 -SNP_pval 1e-6 -minMapQ 30 -minQ 20  -doGlf 2 -out input -nThreads 8
```
- `-doGlf 2` = beagle format
- `-bam BamFileList` = analyze data from bam files 
- `-GL 1`: SAMtools model: GATK assumes sequencing errors are independent, while Samtools believes the second error (unrelated sequencing error that occurs independently of the first one and affects the same genomic position) comes at a higher chance, especially at high depths.
- `-doMajorMinor 1` = infer the major and minor alleles
- `-doMAF 2` = estimate the allele frequencies assuming known minor 
- `-SNP_pval 1e-6` = only keep those sites that have a p-value less than 1e-6 for being variable
  

&nbsp;
&nbsp;

#### 2. View population information file (made previously)
2.1 view a summary of the population information file, cut the first column, sort and count:
```
cut -f 1 -d " " PopInfo | sort | uniq -c
```
2.2 make a population label file
```
cut -f1 -d" " PopInfo > poplabel
```
&nbsp;
&nbsp;
#### 3. View the genotype likelihood beagle file

3.1 View the first 10 columns and 10 lines of the input file:
```
gunzip -c input.beagle.gz | head -n 10 | cut -f 1-10 | column -t
```
3.2 count the number of lines of the input file. The number of lines, indicates the number of loci for which there are GLs plus one (as the command includes the count of the header line):
```
gunzip -c input.beagle.gz | wc -l
```
&nbsp;
&nbsp;


#### 4. Run an analysis of the GLs with NGSadmix, assuming the number of ancestral populations is K:
```
for k in {1..10}; do angsd/NGSadmix -likes input.beagle.gz -K $k -seed 1 -P 25 -o NGSadmix/NGSadmix_$k; done

```
- Trying different values of `-k`, the number of assumed ancestral populations: 2-10
- `-mTol`= include high quality genotypes only
- `-seed 1`= seed for initial guess in EM algorithm
- `-P`= number of threads
- Note: I ran this analysis 3 times to test convergence using `-seed` 1-3.

&nbsp;
&nbsp;

#### 5. Check output file to look at first 5 admixture proportions 

```
head -n 5 NGSadmix_k.qopt
```
&nbsp;
&nbsp;


## Plot the results in R
&nbsp;
#### 6. Data & libraries 
6.1 Libraries
```
library(ggplot2) # Just in case
library(viridis) # Color palette 
```
6.2 Data
```
q2<-read.table("NGSadmix_2.qopt")
pop<-read.csv("PopLabel.csv", header = FALSE)
```
&nbsp;
&nbsp;
#### 7. ID and population info for each individual
```
pop<-as.matrix(pop)
ord = order(pop) #order by population
```
&nbsp;
&nbsp;
#### 8. Plot
8.1 set population parameters and colors
```
num_pops <- 2  # Adjust based on K
colors <- viridis(num_pops)
```
8.2 create structure plot and save as a .png file
```
png(file="K2.png", width = 800, height = 600)
barplot(t(q2)[, ord], col = my_colors, names = pop[ord], las = 2, ylab = "K=2, Admixture Proportions", cex.names = 0.75)
dev.off()
```
&nbsp;


$${\color{lightgreen}Repeat \space steps \space 6-8 \space to \space create \space a \space structure \space plot \space for \space each \space value \space of \space k}$$

&nbsp;
&nbsp;

#### 9. Get an accurate estimate of the “best” K 
9.1 Libraries
```
library(stringr)
```

9.2 Data
```
data<-list.files("structure_outs/", pattern = ".log", full.names = T)

#sanity check
data
```
9.3 use `lapply` to read in all log files at once
```
bigData<-lapply(1:21, FUN = function(i) readLines(data[i]))
```
9.4 pull out the line that starts with "b" from each file and return it as a list
```
foundset<-sapply(1:21, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])

# sanity check
foundset
```
9.5 pull out the first number in the string using the function `sub`
```
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )
```
9.6 store it in a dataframe (index 1:10, corresponding to K values)
```
logs<-data.frame(K = rep(1:7, each=3))
```
9.7 add to likelihood values
```
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))
```
9.8 calculate delta K & probability (use these values to select K, which will be the one that has the highest value)
```
tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))
```


&nbsp;

## EvalAdmix 
&nbsp;
#### 10. Evaluate the fit of admixture models based on estimating the correlation of the residual difference between the true genotypes and the genotypes predicted by the model

```
 ../../../angsd/evalAdmix/evalAdmix -beagle ../BeagleMadeBeforeRagTag/noCLRA.beagle.gz -fname NGSadmix_1.fopt.gz -qname NGSadmix_1.qopt -P 10
```
- *Repeat for K=1-10*


