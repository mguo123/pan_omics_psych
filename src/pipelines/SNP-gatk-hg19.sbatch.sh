#!/bin/bash
#SNP-gatk-hg19.sbatch
#SBATCH --job-name=sample
#SBATCH --output=/oak/stanford/groups/howchang/users/yangzhao/Chang_lab/Diana/human_autoImmune/logs/sample_gatk.out
#SBATCH --error=/oak/stanford/groups/howchang/users/yangzhao/Chang_lab/Diana/human_autoImmune/logs/sample_gatk.err
#SBATCH --time=24:00:00
#SBATCH -p rbaltman,khavari,normal,owners
#SBATCH --nodes=1
#SBATCH --mem=100000
#SBATCH -c 16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mguo1238@stanford.edu
source ~/.bashrc
module load biology
ml gatk
# ml bwa
module load bowtie2
module load samtools
module load bedtools
module load py-macs2
module load tbb
module load java
module load python
module load perl
module load R

BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19/ucsc.hg19.fasta"
# BWA_fasta_ref="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg38/Homo_sapiens_assembly38.fasta"
# BWA_fasta_ref="/oak/stanford/groups/khavari/users/yzhao8/Reference/BWA/hg38_gatk"
 /oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19
data_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bam_files"
output_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp"
refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/hg19"
# # make a merged atac bed file of all the tissue atac peaks
# /oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files]$ cat H9iN-day0.narrowPeak.bed H9iN-day28.narrowPeak.bed H9-Ngn2.narrowPeak.bed SLC_baseline.narrowPeak.bed SLC-Ngn2.narrowPeak.bed SL-Ngn2.narrowPeak.bed H9_baseline.narrowPeak.bed H9iN-day10.narrowPeak.bed H9iN-day4.narrowPeak.bed SL_baseline.narrowPeak.bed  SLC.narrowPeak.bed SL.narrowPeak.bed | bedtools sort -i stdin| bedtools merge -i stdin > all_neuro_pks.bed
interval_file="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files/all_neuro_pks.bed"


#Astrocytes H9_D0 H9_D10 H9_D28 H9_D2 H9_D4 SLC_D0 SLC_D2 SL_D0 SL_D2 
sample=''
# picard_path="/oak/stanford/groups/rbaltman/mguo123/picard/build/libs"

# # list fastq files
# cd $data_path
# fastqList_1=$(find "$(pwd)" -name "sample*_R1_001.fastq.gz" )
# fastqList_2=$(find "$(pwd)" -name "sample*_R2_001.fastq.gz" )
cd $output_path

# # make bam files and index (deduped)
# bwa mem -t 16 -R '@RG\tID:scleroderma\tSM:sample\tPL:illumina\tLB:sample\tPU:sample' $BWA_fasta_ref/Homo_sapiens_assembly38 $fastqList_1 $fastqList_2 > $output_path/sample.sam
# java -jar $picard_path/picard.jar SortSam INPUT=$output_path/sample.sam OUTPUT=$output_path/sample_sorted.bam SORT_ORDER=coordinate
# java -jar $picard_path/picard.jar MarkDuplicates \
#       I=$output_path/sample_sorted.bam \
#       O=$output_path/sample_sorted_mkdup.bam \
#       M=$output_path/sample_marked_dup_metrics.txt
# java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/sample_sorted_mkdup.bam



#First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
gatk BaseRecalibrator \
-nct 80 \
-R $BWA_fasta_ref_file \
-I $data_path/${sample}.bam \
-o $output_path/${sample}_recal.grp \
-knownSites $refSNP_path/dbsnp_138.hg19.vcf \
-knownSites $refSNP_path/1000G_phase1.snps.high_confidence.hg38.vcf \
-knownSites $refSNP_path/1000G_phase1.indels.hg19.sites.vcf \
-knownSites $refSNP_path/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-L $interval_file \
-cov ContextCovariate \
-cov CycleCovariate

# get read haplotype information
gatk PrintReads \
-nct 80 \
-R $BWA_fasta_ref_file \
-I $data_path/${sample}.bam \
-BQSR $output_path/${sample}_recal.grp \
-L $interval_file \
-o $output_path/${sample}_recali.bam \
gatk HaplotypeCaller  \
-nct 80 \
-R $BWA_fasta_ref_file \
-I $output_path/${sample}_recali.bam \
-O $output_path/${sample}_gvcf.gz \
-L $interval_file \
-ERC GVCF








