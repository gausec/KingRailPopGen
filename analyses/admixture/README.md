## Estimating admixture using an ABBABABA test (D-statistic) 

--- 

### Make bam file list
```
ls BAM/*.bam > bam_list.txt
```
### Calculate the allele frequencies
```
angsd -bam bam_list.txt -doMajorMinor 1 -doMaf 1 -GL 2 -out allele.out
```
### ABBABABA test
```
angsd -bam bam_list.txt -doAbbababa 1 -anc CLRA.fasta -max_depth 1000 -min_maf 0.05 -block_size 5000000 -remove_bads 1 -out output
```
- `bam`: Specifies the input bam file list
- `doAbbababa`: Enables the abbababa test
- `anc`: Specifies the reference genome file
- `out`: Specifies the name of output 

