#!/bin/bash

CHR=$1
COUNT=$2

# Hardcoded stuff:
BPSPACE=250000
BASEDIR="/home/juan/tesina/1000Genomes_data"

if [ -z $1 ] || [ -z $2 ]; then
    echo ""
    echo "Please pass two positional arguments:"
    echo "$0 <chromosome_number> <number_of_SNPs_to_keep>"
    echo ""
    exit
fi

OUTDIR="thinned-chromosomes"
OUTFILE="thinned.chr${CHR}.total${COUNT}snps.every${BPSPACE}"

mkdir -p "${BASEDIR}/${OUTDIR}"

plink --vcf ${BASEDIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snps-only --thin-count ${COUNT} --bp-space ${BPSPACE} --recode A-transpose --out ${BASEDIR}/${OUTDIR}/${OUTFILE}

# Keep the 2nd field (dbsnp ID) and then from the 8th until the end (genotypes)
cut ${BASEDIR}/${OUTDIR}/${OUTFILE}.traw -f 2,8-

