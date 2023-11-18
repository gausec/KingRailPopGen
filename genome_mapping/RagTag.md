### Improving reference genome assembly using [RagTag](https://github.com/malonge/RagTag).

- Install
```
conda install -c bioconda ragtag
```
- scaffold utility 
```
ragtag.py scaffold Chicken.fna Rallus_crepitans_1.0.fasta -o CLRA.out
```
- rename output fasta in the CLRA.out directory
```
mv ragtag.scaffold.fasta ../../Improved.CLRA.fasta
```
- extract first 5 chromosomes (based on [chicken chromosome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/))
```
awk -F'\t' 'BEGIN {OFS="\t"} /^CM028482.1_RagTag/ {print $1, $2, $3}' ragtag.scaffold.agp >> Chr1-5.sites.txt; awk -F'\t' 'BEGIN {OFS="\t"} /^CM028483.1_RagTag/ {print $1, $2, $3}' ragtag.scaffold.agp >> Chr1-5.sites.txt; awk -F'\t' 'BEGIN {OFS="\t"} /^CM028484.1_RagTag/ {print $1, $2, $3}' ragtag.scaffold.agp >> Chr1-5.sites.txt; awk -F'\t' 'BEGIN {OFS="\t"} /^CM028485.1_RagTag/ {print $1, $2, $3}' ragtag.scaffold.agp >> Chr1-5.sites.txt; awk -F'\t' 'BEGIN {OFS="\t"} /^CM028486.1_RagTag/ {print $1, $2, $3}' ragtag.scaffold.agp >> Chr1-5.sites.txt
```
