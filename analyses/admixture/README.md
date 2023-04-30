## Estimating admixture using an ABBABABA test (D-statistic) 
---
### Create a file that has all bam files in it
```
angsd -out tryagain -doAbbababa 0 -bam BAM -blockSize 5000000 -anc CLRA.fasta
```
```
angsd -bam bam.txt -doCounts 1 -doAbbababa 1 -anc CLRA.fasta -out output
```
- `bam`: Specifies the input bam file list. You need to replace bamfile.list with the name of your file containing the list of bam files you want to use in the analysis.
- `doAbbababa`: Enables the abbababa analysis. The 1 value indicates that the analysis should be performed.
- `anc`: Specifies the reference genome file in fasta format. You need to replace reference.fasta with the name of your reference genome file.
- `out`: Specifies the prefix of the output files. You need to replace output_prefix with the desired prefix for your output files.
