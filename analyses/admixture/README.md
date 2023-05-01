## Estimating admixture using an ABBABABA test (D-statistic) 

--- 

### Make bam file list
```
ls BAM/*.bam > bam_list.txt
```
### Distribution of genetic variation in a population
```
angsd -bam bam_list.txt -doMajorMinor 1 -doMaf 1 -GL 2 -out allele.out
```
- `-bam`: Specifies the input bam file list
- `-doMajorMinor`: Major and minor alleles based on genotype likelihoods
- `-doMaf`: Site frequency spectrum (SFS) based on major and minor allele frequencies
- `-GL`: Sets genotype likelihood model (2 = GATK).
- `-out`: Output file

### ABBABABA test
```
angsd -bam bam_list.txt -doAbbababa 1 -anc CLRA.fasta -max_depth 1000 -min_maf 0.05 -block_size 5000000 -remove_bads 1 -out output
```
- `-doAbbababa`: Enables the abbababa test
- `-anc`: Specifies the reference genome file

