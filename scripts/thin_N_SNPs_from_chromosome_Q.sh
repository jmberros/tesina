#!/bin/bash

# Hardcoded stuff:
BASEDIR="/home/juan/tesina/1000Genomes_data"
COUNTFILE="/home/juan/tesina/notebook/data/chr_SNP_count_in_galanter"
FACTOR=${1:-1}
BPSPACE=${2:-250000}
OUTDIR=${3:-"thinned-chromosomes"}

mkdir -p "${BASEDIR}/${OUTDIR}"

for CHR_COUNT in `cat "${COUNTFILE}"`; do
    CHR_COUNT=(${CHR_COUNT//,/ })
    CHR=${CHR_COUNT[0]}
    COUNT=$(( ${CHR_COUNT[1]} * $FACTOR ))

    OUTFILE="thinned.chr${CHR}.total${COUNT}snps.every${BPSPACE}"

    plink --vcf ${BASEDIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snps-only --thin-count ${COUNT} --bp-space ${BPSPACE} --make-bed --out ${BASEDIR}/${OUTDIR}/${OUTFILE}

done

# Keep the 2nd field (dbsnp ID) and then from the 8th until the end (genotypes)
# cut ${BASEDIR}/${OUTDIR}/${OUTFILE}.traw -f 2,8-
