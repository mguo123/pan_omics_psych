# Preprocessing Pipeline for ATAC data

1. pull bigwig from kundaje pipeline


Ngn2 data:
/oak/stanford/groups/khavari/users/mguo123/psych_project/cheen_data_ngn2/atac_data/kundaje_output

Baseline data

/oak/stanford/groups/khavari/users/mguo123/cheen_data_baseline/kundaje_output]


New data
/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/kundaje_output




in '../data/processed/washu/atac/'

scp mguo123@login.sherlock.stanford.edu:/oak/stanford/groups/khavari/users/mguo123/psych_project/cheen_data_ngn2/atac_data/kundaje_output/*/call-macs2/shard-0/execution/*.pval.signal.bigwig .

scp mguo123@login.sherlock.stanford.edu:/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/kundaje_output/*/call-macs2_pooled/execution/*.pval.signal.bigwig .



what is copied is

H9-Ngn2_R1.trim.merged.nodup.tn5.pval.signal.bigwig
RTTA_R1.trim.merged.nodup.tn5.pval.signal.bigwig
SLC-Ngn2_R1.trim.merged.nodup.tn5.pval.signal.bigwig
SL-Ngn2_R1.trim.merged.nodup.tn5.pval.signal.bigwig
H9iN-day0-A-GCTACGCT_S3_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig
H9iN-day10-A-TAAGGCGA_S5_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig
H9iN-day28-A-CGTACTAG_S6_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig
H9iN-day4-A-AAGAGGCA_S4_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig
SL-A-GGACTCCT_S1_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig
SLC-A-CTCTCTAC_S2_L004_R1_001.trim.merged.nodup.tn5.pooled.pval.signal.bigwig



-------------



if need to create Washu Bigwig tracks
on sherlock







#ml biology
#ml samtools
#ml py-deeptools


#cd /oak/stanford/groups/khavari/users/mguo123/psych_project

#for f in `ls new_atac_021820/merged_bam/*_merged.bam`; do echo $f; filename=${f##*/};sample=${filename%%_*}; echo $sample; bamCoverage -b $f -o atac_visualization/${sample}.bigWig; done


#for f in `ls cheen_data_ngn2/atac_data/atac_footprinting/merged_bam/*bam`; do echo $f; filename=${f##*/};sample=${filename%%_*}; echo $sample; bamCoverage -b $f -o atac_visualization/${sample}.bigWig; done


#cheen_data_ngn2/atac_data/atac_footprinting/merged_bam/*bam
