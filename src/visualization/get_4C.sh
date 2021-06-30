#!/bin/sh

# LOCAL VERSION
# get_4c.sh
# enter in gene symbol chrom and position and runs through all samples in this folder to get 4C matrix

cd /oak/stanford/groups/khavari/users/mguo123/psych_project/data/hicpro_matrix_abs_files

samples=(Astro_B1_5000 Astro_B2_5000 H9-B1_5000 H9-B2_5000 H9D10-B1_5000 H9D10-B2_5000 H9D28-B1_5000 H9D28-B2_5000 SL-B1_5000 SL-B2_5000 SLC-B1_5000 SLC-B2_5000)


while getopts "c:p:s:" flag
do
    case "${flag}" in
        c) chr=${OPTARG};;
        p) position=${OPTARG};;
        s) sym=${OPTARG};;
    esac
done
echo "chr: $chr";
echo "position: $position";
echo "sym: ${sym}";
#chr=4
#position=113435000
#sym=NEUROG2



for sample in "${samples[@]}"; do
echo $sample

mat_file=$sample.matrix.bed

#chr=4
#position=113435000
#sym=NEUROG2



awk -v position=$position -v chr=$chr  'BEGIN{OFS="\t"} {if ($1==chr && ($2==position|| $5==position)) {if ($2==position){print $1,$2,$3,$4,$5,$6,$7} else {print $4,$5,$6,$1,$2,$3,$7}};}' $mat_file >  4C_files/4C_${sym}_${sample}.matrix.bed



done
