### Improving reference genome assembly using [RagTag](https://github.com/malonge/RagTag).

- Install
```
conda install -c bioconda ragtag
```
- scaffold utility 
```
ragtag.py scaffold Chicken.fna Rallus_crepitans_1.0.fasta -o CLRA.out
```
-rename output fasta in the CLRA.out directory
```
mv ragtag.scaffold.fasta ../../Improved.CLRA.fasta
```
