### I will be using the [DownsampleSam tool](https://broadinstitute.github.io/picard/) within [Picard tools](https://broadinstitute.github.io/picard/) from [GATK](https://gatk.broadinstitute.org/hc/en-us) to downsample bam files for downstream analyses.
&nbsp; 
&nbsp; 


#### 1. Set up

&nbsp; 1.1 Get newest verion of java
```
wget https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz
```
&nbsp; 1.2 extract
```
tar xvf jdk-21_linux-x64_bin.tar.gz
```
&nbsp; 1.3 shortcut to newest verion of java
```
java=~/jdk-21.0.1/bin
```
&nbsp; 1.4 shortcut to bam files
```
bam=~/WholeGenomes/SortedBAM/
```
&nbsp; 1.5 Get newest version of picard tools (I but this in the java bin directory)
```
wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

```
&nbsp;

#### 2. Downsample all bam files keeping 50% of all reads
```
 java -jar jdk-21.0.1/bin/picard.jar DownsampleSam -I WholeGenomes/SortedBAM/11101.sorted.bam -O 11101.downsampled.bam -P 0.5
```

&nbsp;

&nbsp;

---
&nbsp;


### I am also interested in filtering out chromosomes & IDing others
&nbsp;


#### 1. Make blast databases
1.1 W chromosome
```
makeblastdb -in  W-chr.fasta  -dbtype nucl -out W_chromosome_db
```
1.2 Z chromosome
```
makeblastdb -in Z-chr.fasta -dbtype nucl -out Z_chromosome_db
```
1.3 chromosome 1
```
makeblastdb -in  Chromosome1.fasta  -dbtype nucl -out chr1_db
```
1.4 chromosome 2
```
makeblastdb -in  Chromosome2.fasta  -dbtype nucl -out chr2_db
```
1.5 chromosome 3
```
makeblastdb -in  Chromosome3.fasta  -dbtype nucl -out chr3_db
```
1.6 chromosome 4
```
makeblastdb -in  Chromosome4.fasta  -dbtype nucl -out chr4_db
```
1.7 chromosome 5
```
makeblastdb -in  Chromosome5.fasta  -dbtype nucl -out chr5_db
```
&nbsp;
#### 2. Perform Sequence Searches against reference 
1.1 W chromsome
```
blastn -query ../11101.fasta -db W_chromosome_db -out blast_output/W_chromosome_results.tx
```
1.2 Z chromosome
```
blastn -query ../11101.fasta -db Z_chromosome_db -out blast_output/Z_chromosome_results.txt
```
1.2 chromosome 1
```
blastn -query ../11101.fasta -db chr1_db -out -evalue 1e-6 blast_output/chromosome1_results.txt
```
1.3 chromosome 2
```
blastn -query ../11101.fasta -db chr2_db -out -evalue 1e-6 blast_output/chromosome2_results.txt
```
1.4 chromosome 3
```
blastn -query ../11101.fasta -db chr3_db -out -evalue 1e-6 blast_output/chromosome3_results.txt
```
1.5 chromosome 4
```
blastn -query ../11101.fasta -db chr4_db -out -evalue 1e-6 blast_output/chromosome4_results.txt
```
1.6 chromosome 5
```
blastn -query ../11101.fasta -db chr5_db -out -evalue 1e-6 blast_output/chromosome5_results.txt
```


