### I will be using [Picard tools](https://broadinstitute.github.io/picard/) form GATK to downsample bam files for downstream analyses.
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

#### 2. Downsample all bam files keeping 10% of all reads
```
$java java -v jdk-21.0.1 -jar picard.jar DownsampleSam \
       I=input.bam \
       O=downsampled.bam \
       P=0.1
