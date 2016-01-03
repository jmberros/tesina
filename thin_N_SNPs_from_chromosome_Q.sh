#!/bin/bash

# Hardcoded stuff:
BPSPACE=250000
BASEDIR="/home/juan/tesina/1000Genomes_data"
COUNTFILE="/home/juan/tesina/notebook/data/chr_SNP_count_in_galanter"
OUTDIR="thinned-chromosomes"

mkdir -p "${BASEDIR}/${OUTDIR}"

for CHR_COUNT in `cat "${COUNTFILE}"`; do
    CHR_COUNT=(${CHR_COUNT//,/ })
    CHR=${CHR_COUNT[0]}
    COUNT=${CHR_COUNT[1]}

    OUTFILE="thinned.chr${CHR}.total${COUNT}snps.every${BPSPACE}"

    plink --vcf ${BASEDIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snps-only --thin-count ${COUNT} --bp-space ${BPSPACE} --make-bed --out ${BASEDIR}/${OUTDIR}/${OUTFILE}

done

# Keep the 2nd field (dbsnp ID) and then from the 8th until the end (genotypes)
# cut ${BASEDIR}/${OUTDIR}/${OUTFILE}.traw -f 2,8-
