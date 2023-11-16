### I am interested in filtering out chromosomes & IDing chromosomes to keep
---
&nbsp;


#### 1. Make blast databases
1.1 W chromosome
```
makeblastdb -in  W-chr.fasta  -dbtype nucl -out W_chromosome_db
```
&nbsp;

1.2 Z chromosome
```
makeblastdb -in Z-chr.fasta -dbtype nucl -out Z_chromosome_db
```
&nbsp;

1.3 chromosome 1
```
makeblastdb -in  Chromosome1.fasta  -dbtype nucl -out chr1_db
```
&nbsp;

1.4 chromosome 2
```
makeblastdb -in  Chromosome2.fasta  -dbtype nucl -out chr2_db
```
&nbsp;

1.5 chromosome 3
```
makeblastdb -in  Chromosome3.fasta  -dbtype nucl -out chr3_db
```
&nbsp;

1.6 chromosome 4
```
makeblastdb -in  Chromosome4.fasta  -dbtype nucl -out chr4_db
```
&nbsp;

1.7 chromosome 5
```
makeblastdb -in  Chromosome5.fasta  -dbtype nucl -out chr5_db
```
&nbsp;
#### 2. BAM query file > FASTA (Note: this wasn't the correct format, so I used the reference fasta as a query)
```
samtools fasta 11101.downsampled.bam > 11101.fasta 
```

&nbsp;

#### 3. Perform Sequence Searches against reference 
3.1 W chromsome
```
blastn -query 11101.fasta -db W_chromosome_db -evalue 1e-2 -out blast_output/W_chromosome_results.tx
```
&nbsp;

3.2 Z chromosome
```
blastn -query 11101.fasta -db Z_chromosome_db -evalue 1e-2 -out blast_output/Z_chromosome_results.txt
```
&nbsp;

3.3 chromosome 1
```
blastn -query Rallus_crepitans_1.0.fasta -db chr1_db -evalue 1e-6 -out blast_output/chromosome1_results.txt
```
&nbsp;

3.4 chromosome 2
```
blastn -query Rallus_crepitans_1.0.fasta -db chr2_db -evalue 1e-6 -out blast_output/chromosome2_results.txt
```
&nbsp;

3.5 chromosome 3
```
blastn -query Rallus_crepitans_1.0.fasta -db chr3_db -evalue 1e-6 -out blast_output/chromosome3_results.txt
```
&nbsp;

3.6 chromosome 4
```
blastn -query Rallus_crepitans_1.0.fasta -db chr4_db -evalue 1e-6 -out blast_output/chromosome4_results.txt
```
&nbsp;

3.7 chromosome 5
```
blastn -query Rallus_crepitans_1.0.fasta -db chr5_db -evalue 1e-6 -out blast_output/chromosome5_results.txt
```
&nbsp;
&nbsp;


#### 4. Extract HSP info
```
perl extract_HSP_information.pl chromosome1_results.txt
```
*Repeat for all output files*
&nbsp;

*This script can be found in the bin subdirectory*

