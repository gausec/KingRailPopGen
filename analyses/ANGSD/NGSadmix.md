## NGSadmix
&nbsp;
> I am using [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmixTutorial) in ANGSD to estimate individual admixture proportions from NGS data for  four populations of king rails, and then examine whether these populations separate based on their admixture patterns. I am looking for substructure among the populations.

&nbsp;
&nbsp;
#### 1. Create the the beagle genotype likelihood input file using ANGSD: calculate genotype likelihoods for polymorphic sites using ANGSD (NGSadmix uses Genotype Likelihoods (GLs) in .beagle format as input)
1.1 shortcut
```
ANGSD=../../angsd
```
1.2 Beagle file
```
$ANGSD/angsd -bam BamFileList.txt -ref CLRAindex/Rallus_crepitans_1.0.fasta -anc CLRAindex/Rallus_crepitans_1.0.fasta -GL 1 -doSaf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -minind 5 -SNP_pval 1e-6 -minMapQ 30 -minQ 20  -doGlf 2 -out input.gz -nThreads 8
```
- `-GL 1` = SAMtools model
- `-doGlf 3` = beagle binary
- `-doGlf 2` = beagle format

&nbsp;
&nbsp;

#### 2. View population information file (made previously)
2.1 view a summary of the population information file, cut the first column, sort and count:
```
cut -f 1 -d " " PopInfo | sort | uniq -c
```
make a population label file and place it in the output directory
```
cut -f1 -d" " PopInfo > poplabel
```
&nbsp;
&nbsp;
#### 3. View the genotype likelihood beagle file

3.1 View the first 10 columns and 10 lines of the input file:
```
gunzip -c input.gz.beagle.gz | head -n 10 | cut -f 1-10 | column -t
```
3.2 count the number of lines of the input file. The number of lines, indicates the number of loci for which there are GLs plus one (as the command includes the count of the header line):
```
gunzip -c input.gz.beagle.gz | wc -l
```
&nbsp;
&nbsp;


#### 4. Run an analysis of the GLs with NGSadmix, assuming the number of ancestral populations is K:
```
 for ((k=1; k<=10; k++)); do $ANGSD/NGSadmix -likes input.gz.beagle.gz -K $k -minMaf 0.05 -seed 1 -P 10 -tolLike50 -o NGSadmix_$k; done
```
- Trying `-k` = 1-10


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
q<-read.table("NGSadmix_k.qopt")
pop<-read.csv("PopLabel.csv", header = FALSE)
```
&nbsp;
&nbsp;
#### 7. ID and population info for each individual
```
pop<-as.matrix(pop)
ord = order(pop)
```
&nbsp;
&nbsp;
#### 8. Plot
```
num_pops <- 1  # Adjust based on K
colors <- viridis(num_pops)
barplot(t(q)[, ord], col = colors, names = pop[ord], las = 2, ylab = "Demo1 Admixture Proportions", cex.names = 0.75)
```
&nbsp;
&nbsp;
