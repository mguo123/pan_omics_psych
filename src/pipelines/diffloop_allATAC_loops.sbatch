#!/bin/bash
#
#SBATCH --job-name=diffLoops_allATAC_Loops_psych
#SBATCH --output=/oak/stanford/groups/khavari/users/yzhao8/Margaret/GWAS/logs/diffLoops_allATAC_Loops_psych.out
#SBATCH --error=/oak/stanford/groups/khavari/users/yzhao8/Margaret/GWAS/logs/diffLoops_allATAC_Loops_psych.err
#SBATCH --time=48:00:00
#SBATCH --partition=khavari,rbaltman,normal,owners
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mguo123@stanford.edu

cd /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment
ml biology bedtools

###### STEP 1: make diff atac files from diffloop bedpe and atac bed files
# (base) [mguo123@sh02-10n18 /oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files]$ ls
# Astrocytes_merged.bed  H9_D0_merged.bed  H9_D10_merged.bed  H9_D28_merged.bed  H9_D2_merged.bed  H9_D4_merged.bed  old_files  SLC_D0_merged.bed  SLC_D2_merged.bed  SL_D0_merged.bed  SL_D2_merged.bed


# (base) [mguo123@sh02-10n18 /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment/diffloop_allatac_anchors]$ ls  ../../hichip_diffloop_files/*anchor.bed
# ../../hichip_diffloop_files/Astro_diffloop_anchor.bed  ../../hichip_diffloop_files/H9_D10_diffloop_anchor.bed  ../../hichip_diffloop_files/H9_D2_diffloop_anchor.bed   ../../hichip_diffloop_files/SLC_D2_diffloop_anchor.bed  ../../hichip_diffloop_files/SL_D2_diffloop_anchor.bed
# ../../hichip_diffloop_files/H9_D0_diffloop_anchor.bed  ../../hichip_diffloop_files/H9_D28_diffloop_anchor.bed  ../../hichip_diffloop_files/SLC_D0_diffloop_anchor.bed  ../../hichip_diffloop_files/SL_D0_diffloop_anchor.bed

# tissue="Astro"
# tissue="H9_D0" 
# tissue="H9_D10" 
# tissue="H9_D28" 
# tissue="H9_D2" 
# tissue="SL_D0" 
# tissue="SL_D2" 
# tissue="SLC_D0" 
# tissue="SLC_D2" 

# tissue="SL_D2" 
# atac_tissue=${tissue}
# #except Astrocytes!!!
# cd /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment/diffloop_allatac_anchors
# awk '{print "chr"$1"\t"$2"\t"$3}' ../../hichip_diffloop_files/${tissue}_diffloop_anchor.bed | sort  -k1,1 -k2,2n| uniq |  bedtools intersect -a stdin -b ../../atac_bed_files/${atac_tissue}_merged.bed -wa | uniq > ${tissue}_diffloop_atac_anchor.bed 
# cd ../diffloop_allatac_loop
# awk '{print "chr"$1"\t"$2"\t"$3"\t""chr"$4"\t"$5"\t"$6"\t"$7}' ../../hichip_diffloop_files/${tissue}_diffloop_loop.bed | bedtools intersect -a stdin -b ../diffloop_allatac_anchors/${tissue}_diffloop_atac_anchor.bed  -wa > ${tissue}_diffloop_atac_loop.bed 


# (base) [mguo123@sh02-10n18 /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment]$ wc -l diffloop_allatac_anchors/*
#   13952 diffloop_allatac_anchors/Astro_diffloop_atac_anchor.bed
#   32687 diffloop_allatac_anchors/H9_D0_diffloop_atac_anchor.bed
#   19653 diffloop_allatac_anchors/H9_D10_diffloop_atac_anchor.bed
#   13694 diffloop_allatac_anchors/H9_D28_diffloop_atac_anchor.bed
#       0 diffloop_allatac_anchors/H9_D2_diffloop_atac_anchor.bed
#   40071 diffloop_allatac_anchors/SLC_D0_diffloop_atac_anchor.bed
#       0 diffloop_allatac_anchors/SLC_D2_diffloop_atac_anchor.bed
#   38638 diffloop_allatac_anchors/SL_D0_diffloop_atac_anchor.bed
#       0 diffloop_allatac_anchors/SL_D2_diffloop_atac_anchor.bed
#  158695 total

## NOTE: SKIPP INg the NGN2s

# ####### STEP 2: make blacklist exclusion files
# cd /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment
# tissues=( Astro  H9_D0 H9_D10 H9_D28 H9_D2 SL_D0 SL_D2 SLC_D0 SLC_D2)
# for tissue in "${tissues[@]}";
# do
# echo $tissue
# cat ./resources/hg19.blackList.bed ./diffloop_allatac_anchors/${tissue}_diffloop_atac_anchor.bed  > ./diffloop_atac_anchors_exclude/${tissue}_excl.bed
# done


###### STEP 3: shuffling
#results dir
mkdir -p result_diffloop_atac

cd /oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment
tissues=( Astro  H9_D0 H9_D10 H9_D28 SL_D0  SLC_D0)
for tissue in "${tissues[@]}";
do
echo $tissue

#REAL THING
 cut -f1-3 ./diffloop_allatac_loop/${tissue}_diffloop_atac_loop.bed > ./result_diffloop_atac/${tissue}_ACTUAL_anchor1.bed
 cut -f4-6 ./diffloop_allatac_loop/${tissue}_diffloop_atac_loop.bed > ./result_diffloop_atac/${tissue}_ACTUAL_anchor2.bed
 cat ./result_diffloop_atac/${tissue}_ACTUAL_anchor1.bed ./result_diffloop_atac/${tissue}_ACTUAL_anchor2.bed > ./result_diffloop_atac/${tissue}_ACTUAL_anchors.bed
 rm ./result_diffloop_atac/${tissue}_ACTUAL_anchor1.bed
 rm ./result_diffloop_atac/${tissue}_ACTUAL_anchor2.bed
bedtools intersect -a ./resources/LinkedSNP-all-sorted.bed -b ./diffloop_allatac_anchors/${tissue}_diffloop_atac_anchor.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_ACTUAL_overlap_all.bed
bedtools intersect -a ./resources/psych_mpra_lib.bed -b ./diffloop_allatac_anchors/${tissue}_diffloop_atac_anchor.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_ACTUAL_overlap_psychmpra.bed
bedtools intersect -a ./resources/SNP-nonCancer-nonFilter-1.bed  -b ./diffloop_allatac_anchors/${tissue}_diffloop_atac_anchor.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_ACTUAL_overlap_noncancer.bed
 rm ./result_diffloop_atac/${tissue}_ACTUAL_anchors.bed

#SHUFFLED THING
for m in {1..100}
do
 bedtools shuffle -i ./diffloop_allatac_loop/${tissue}_diffloop_atac_loop.bed -g ./resources/chrom_hg19.sizes -seed $m -excl ./diffloop_atac_anchors_exclude/${tissue}_excl.bed -bedpe > ./result_diffloop_atac/${tissue}_"$m".bed
 cut -f1-3 ./result_diffloop_atac/${tissue}_"$m".bed > ./result_diffloop_atac/${tissue}_"$m"_anchor1.bed
 cut -f4-6 ./result_diffloop_atac/${tissue}_"$m".bed > ./result_diffloop_atac/${tissue}_"$m"_anchor2.bed
 cat ./result_diffloop_atac/${tissue}_"$m"_anchor1.bed ./result_diffloop_atac/${tissue}_"$m"_anchor2.bed > ./result_diffloop_atac/${tissue}_"$m"_anchors.bed
 rm ./result_diffloop_atac/${tissue}_"$m"_anchor1.bed
 rm ./result_diffloop_atac/${tissue}_"$m"_anchor2.bed
 bedtools intersect -a ./resources/LinkedSNP-all-sorted.bed -b ./result_diffloop_atac/${tissue}_"$m"_anchors.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_"$m"_overlap_all.bed
 bedtools intersect -a ./resources/psych_mpra_lib.bed -b ./result_diffloop_atac/${tissue}_"$m"_anchors.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_"$m"_overlap_psychmpra.bed
 bedtools intersect -a ./resources/SNP-nonCancer-nonFilter-1.bed  -b ./result_diffloop_atac/${tissue}_"$m"_anchors.bed -wa | sort -k1,1 -k2,2n | uniq > ./result_diffloop_atac/${tissue}_"$m"_overlap_noncancer.bed
 rm ./result_diffloop_atac/${tissue}_"$m"_anchors.bed
#done shuffle loop
done


# done tissues loop
done
