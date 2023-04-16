### Variant Calling

### Steps
1. Using
```
bcftools mpileup -f Rallus_crepitans_1.0.fasta -Ou -q 20 -Q 20 -C 50 -I -d 1000000 -a AD $file.sorted.bam | bcftools call -mv -Ov -f GQ -o $file_output_raw.vcf
```
