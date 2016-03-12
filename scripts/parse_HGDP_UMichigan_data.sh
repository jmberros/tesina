#!/bin/bash

DIR=$1

for FILE in `ls $DIR*.geno.gz`; do
    zcat $FILE | head -n 5 > ${FILE:0:-8}.map
    zcat $FILE | tail -n +6 > ${FILE:0:-8}.geno
done

