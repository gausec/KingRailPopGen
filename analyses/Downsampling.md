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


### I am also interested in filtering out the sex chromosomes (Z & W)
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
&nbsp;
#### 2. Perform Sequence Searches against reference genome
1.1 W chromsome
```
blastn -query SortedBAM/CLRAindex/Rallus_crepitans_1.0.fasta -db W_chromosome_db -out W_chromosome_results.tx
```
1.2 Z chromosome
```
blastn -query SortedBAM/CLRAindex/Rallus_crepitans_1.0.fasta -db Z_chromosome_db -out Z_chromosome_results.txt
```

&nbsp;

#### 3. 
