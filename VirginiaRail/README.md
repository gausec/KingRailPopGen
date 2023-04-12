# Aligning VIRA genome to VIRA genome (
1. Two Virginina rail Whole gemones from Mackay Island NWR. They have been preprocessed in FastP (filtered and trimmed).
2. Indexing the Virgina rail whole genome from California (Accession no. JAKCOZ000000000) using BWA to use as a reference genome.
```
bwa index GCA_022605955.1_bRalLim1.0.p_genomic.fna.gz
```   
3. We moved all the index files into a new directory  
```
mkdir VIRAindex
```   
4. Aligning the two genomes to a Virgina rail whole genome from California using BWA.
