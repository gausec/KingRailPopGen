# Phylogenetic comparison of three _Rallus_ species
---
### Aim  
We are using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD#Overview) to phlyogenetically compare three species of North American rails, the clapper rail (_Rallus crepitans_) the Virginia rail (_Rallus limicola_) and the king rail (_Rallus elegans_). One individual from each species will be compared using [ABBABABA(D-stat)](http://www.popgen.dk/angsd/index.php/Abbababa). We are using a Virgina rail whole genome from California (Accession no. JAKCOZ000000000) as a reference to align our two Virigia rail genomes. We will choose one of these two individuals pending coverage information.

---
### Contributors:  Carol Gause and Megan Linke  
---

1. Two Virginina rail whole genomes from Mackay Island NWR have been preprocessed using fastp (filtered and trimmed). These were cleaned with all other sequences (see fastp folder in parent directory).  
```
fastp -i 67438_S1_L002_R1_001.fastq.gz -I 67438_S1_L002_R1_001.fastq.gz -o cleaned/67438-R1.fastq.gz -O cleaned/67438-R2.fastq.gz --html reports/67438.html --json reports/67438.json
```
```
fastp -i 78832_S1_L002_R1_001.fastq.gz -I 78832_S1_L002_R1_001.fastq.gz -o cleaned/78832-R1.fastq.gz -O cleaned/78832-R2.fastq.gz --html reports/78832.html --json reports/78832.json
``` 
2. Re-naming and indexing the reference genome using BWA.   
```
mv GCA_022605955.1_bRalLim1.0.p_genomic.fna.gz VIRAref.fna.gz 
```  
```
bwa index VIRAref.fna.gz
```       
3. We moved all the index files into a new directory  
```
mkdir VIRAindex 
``` 
```
mv VIRAref.fna.gz /VIRAindex #Repeated for all index files
```  
4. Align the two genomes to the Virgina rail reference genome using BWA and convert it to a BAM file using [SAMtools](https://github.com/samtools/samtools). 
```
for r1file in *R1.fastq.gz; do base=$(basename -s -R1.fastq.gz ${r1file}); bwa mem -t 40 VIRAindex/VIRAref.fna.gz ${base}-R1.fastq.gz ${base}-R2.fastq.gz | samtools view -@ 40 -bS > SAM/${base}.bam; done 
```
5. Arrange the aligned reads in the BAM files based on their coordinates (Sorting) using SAMtools.
```
mkdir sorted
mkdir temp_dir
```
```
for file in 67438.bam 78832.bam; do echo "Sorting ${file}..."; samtools sort -@ 40 -T temp_dir ${file} -o sorted/${file%.*}.sorted.bam; done
```
6. Determine the average depth of the sorted bam files
```
samtools depth 67438.sorted.bam | awk '{sum+=$3} END {print "Average depth:", sum/NR}'
```
*The average depth for 67438 is 9.5x*
```
samtools depth 78832.sorted.bam | awk '{sum+=$3} END {print "Average depth:", sum/NR}'
```
*The average depth for 78832 is 6.6x.*
