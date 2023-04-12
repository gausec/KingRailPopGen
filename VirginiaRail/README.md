# Phylogenetic comparison of three _Rallus_ species
---
### Aim
I am using [ANGSD](http://www.popgen.dk/angsd/index.php/Abbababa) to phlyogenetically compare three species of North American Rails, the Clapper rail (_Rallus _) the Virginia rail (_Rallus limicola_) and the king rail (_Rallus elegans_). One individual from each species will be compared using ABBABABA(D-stat).
---

### Contributors: Carol Gause and Megan Linke
---

1. Two Virginina rail Whole gemones from Mackay Island NWR. They have been preprocessed in FastP (filtered and trimmed).
2. Indexing the Virgina rail whole genome from California (Accession no. JAKCOZ000000000) using BWA to use as a reference genome.
```
bwa index GCA_022605955.1_bRalLim1.0.p_genomic.fna.gz
```   
3. We moved all the index files into a new directory  
```
mkdir VIRAindex
```   
4. Align the two genomes to the Virgina rail reference genome using BWA and convert it to a BAM file using [SAMtools](https://github.com/samtools/samtools). 
5. Arrange the aligned reas in the BAM files based on their coordinates (Sorting) using SAMtools.

