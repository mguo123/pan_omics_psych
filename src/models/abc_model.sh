#!/bin/bash
### psych project

# get modules
ml biology samtools bedtools
ml devel java
export PATH=/oak/stanford/groups/khavari/users/mguo123/juicer/SLURM/scripts:$PATH
ml system curl zlib

# go to working directory
cd /oak/stanford/groups/khavari/users/mguo123/psych_project/code/ABC-Enhancer-Gene-Prediction

# establish folders
atac_bed_dir=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files/
atac_bam_dir=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bam_files/
reference_dir=/oak/stanford/groups/khavari/users/mguo123/Reference/
rna_dir=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/rna_tpm_files/
hichip_dir=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/hichip_fithichip_bedpe_files/
tss_filepath=/oak/stanford/groups/khavari/users/mguo123/psych_project/annon/gencode.v19_bed6.gene.bed

# loop through each tissue
# tissue=H9_D0 #other tissues here
for tissue in Astrocytes  H9_D0  H9_D10  H9_D2  H9_D28  SLC_D0  SLC_D2  SL_D0  SL_D2;
do 
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo $tissue
out_dir=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/ABC_out/${tissue}/
mkdir -p $out_dir
mkdir -p ${out_dir}Peaks
mkdir -p ${out_dir}Neighborhoods
mkdir -p ${out_dir}Predictions


# Step 1 - get candidate regions
python src/makeCandidateRegions.py \
--narrowPeak ${atac_bed_dir}${tissue}_merged.bed \
--bam ${atac_bam_dir}${tissue}.bam \
--outDir ${out_dir}Peaks/ \
--chrom_sizes ${reference_dir}RSEM/hg19/hg19.chrom.sizes \
--regions_blacklist ${reference_dir}regions/wgEncodeDacMapabilityConsensusExcludable.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000 

# Step 2 - quantify enhancer activity
python src/run.neighborhoods.py \
--candidate_enhancer_regions ${out_dir}/Peaks/${tissue}_merged.bed.candidateRegions.bed \
--genes ${tss_filepath} \
--ATAC ${atac_bam_dir}${tissue}.bam \
--expression_table ${rna_dir}${tissue}.TPM.txt \
--chrom_sizes ${reference_dir}RSEM/hg19/hg19.chrom.sizes \
--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType ${tissue} \
--outdir ${out_dir}Neighborhoods/ 

# Step 3: Predict
python src/predict.py \
--enhancers ${out_dir}Neighborhoods/EnhancerList.txt \
--genes ${out_dir}Neighborhoods/GeneList.txt \
--HiCdir ${hichip_dir}${tissue}/ \
--chrom_sizes ${reference_dir}RSEM/hg19/hg19.chrom.sizes \
--hic_type bedpe \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType ${tissue} \
--outdir ${out_dir}Predictions/ \
--make_all_putative

done