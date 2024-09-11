#!/bin/bash

# paths
ANGSD_DIR="../../angsd/misc"
REF_FASTA="../../ImprovedCLRAindex/ImprovedCLRA.fasta"

# samples
SAMPLES=("NC" "LA" "FL" "OH")

# Folded SFS
for SAMPLE in "${SAMPLES[@]}"; do
    $ANGSD_DIR/realSFS ${SAMPLE}.saf.idx -P 24 -fold 1 -anc $REF_FASTA -ref $REF_FASTA > ${SAMPLE}.Folded.SFS
done

# Folded SFS to theta estimates
for SAMPLE in "${SAMPLES[@]}"; do
    $ANGSD_DIR/realSFS saf2theta ${SAMPLE}.saf.idx -outname ${SAMPLE} -sfs ${SAMPLE}.Folded.SFS -fold 1
done

# theta statistics
for SAMPLE in "${SAMPLES[@]}"; do
    $ANGSD_DIR/thetaStat do_stat ${SAMPLE}.thetas.idx
done

# extract specific values and compute averages
for SAMPLE in "${SAMPLES[@]}"; do
    awk 'NR > 1 {print $5 / $14}' ${SAMPLE}.thetas.idx.pestPG >> ${SAMPLE}pi.txt
    awk 'NR > 1 {print $4 / $14}' ${SAMPLE}.thetas.idx.pestPG >> ${SAMPLE}tW.txt
done

# final averages
for SAMPLE in "${SAMPLES[@]}"; do
    awk '{ sum += $1 } END { print sum / NR }' ${SAMPLE}pi.txt >> ${SAMPLE}pi.final
    awk '{ sum += $1 } END { print sum / NR }' ${SAMPLE}tW.txt >> ${SAMPLE}tW.final
done
