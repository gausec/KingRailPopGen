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
&nbsp; 1.5 Get newest version of picard tools (in the java bin directory)
```
wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

```
&nbsp;

#### 2. Example code: downsample a given bam file keeping 80% of all reads
```
 java -jar jdk-21.0.1/bin/picard.jar DownsampleSam -I WholeGenomes/SortedBAM/11101.sorted.bam -O 11101.downsampled.bam -P 0.80
```

&nbsp;

&nbsp;

#### 3. Sanity check using [samtools quickcheck](https://www.htslib.org/doc/samtools-quickcheck.html)
```
samtools quickcheck -v *.bam > bad_bams.fofn   && echo 'all files passed the check' || echo 'some files failed the check, see bad_bams.fofn'
```
Notes about sanity check script:
- The `-v` flag enables verbose output, which provides more details about any issues.
- The `.fofn` extension stands for "File of File Names." It's a convention used to indicate that the file contains a list of file names, often one per line. Here, `bad_bams.fofn` will contain a list of any files that failed the samtools quickcheck.
- The double pipe `||` in shell scripting runs the second command only if the first command fails. It’s a way to handle errors or run an alternative action if the first command doesn’t work.
&nbsp;
