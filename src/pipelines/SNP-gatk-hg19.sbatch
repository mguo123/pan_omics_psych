#!/bin/bash
#SNP-gatk-hg19-atac.sbatch
#SBATCH --job-name=sample
#SBATCH --output=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp/sample_gatk.out
#SBATCH --error=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp/sample_gatk.err
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
# module load bowtie2
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
 # /oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19
# data_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bam_files"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
# 
output_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp"
refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/hg19"
# # make a merged atac bed file of all the tissue atac peaks
# /oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files]$ cat H9iN-day0.narrowPeak.bed H9iN-day28.narrowPeak.bed H9-Ngn2.narrowPeak.bed SLC_baseline.narrowPeak.bed SLC-Ngn2.narrowPeak.bed SL-Ngn2.narrowPeak.bed H9_baseline.narrowPeak.bed H9iN-day10.narrowPeak.bed H9iN-day4.narrowPeak.bed SL_baseline.narrowPeak.bed  SLC.narrowPeak.bed SL.narrowPeak.bed | bedtools sort -i stdin| bedtools merge -i stdin > all_neuro_pks.bed
interval_file="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files/all_neuro_pks.bed"


#Astrocytes H9_D0 H9_D10 H9_D28 H9_D2 H9_D4 SLC_D0 SLC_D2 SL_D0 SL_D2 
sample='Astrocytes'
picard_path="/oak/stanford/groups/rbaltman/mguo123/picard/build/libs"

# STEP 1(IF START FROM BAM-SKIP)
# # list fastq files
cd $data_path
fastqList_1=$(find "$(pwd)" -name "sample*_R1_001.fastq.gz" )
fastqList_2=$(find "$(pwd)" -name "sample*_R2_001.fastq.gz" )
cd $output_path

# # STEP 2 BAM PREPROCESSING (TAKES A WHILE) make bam files and index (deduped)
# add read group, dedup (already done), sort index
# bwa mem -t 16 -R '@RG\tID:scleroderma\tSM:sample\tPL:illumina\tLB:sample\tPU:sample' $BWA_fasta_ref/Homo_sapiens_assembly38 $fastqList_1 $fastqList_2 > $output_path/sample.sam
# java -jar $picard_path/picard.jar SortSam INPUT=$output_path/sample.sam OUTPUT=$output_path/sample_sorted.bam SORT_ORDER=coordinate
# java -jar $picard_path/picard.jar MarkDuplicates \
#       I=$output_path/sample_sorted.bam \
#       O=$output_path/sample_sorted_mkdup.bam \
#       M=$output_path/sample_marked_dup_metrics.txt
java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
    I=${data_path}/${sample}.bam  \
    O=$output_path/${sample}.bam \
    SORT_ORDER=coordinate \
    RGID=psych \
    RGLB=sample \
    RGPU=sample \
    RGPL=illumina \
    RGSM=${sample} 
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}.bam
## DEBUG CHECK Readgroups
# samtools view -H $output_path/${sample}.bam  | grep '@RG'


# STEP 3: base recalibration takes 5-10 min
#First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
gatk BaseRecalibrator \
-R $BWA_fasta_ref_file \
-I $output_path/${sample}.bam \
-O $output_path/${sample}_recal.table \
--known-sites $refSNP_path/dbsnp_138.hg19.vcf \
--known-sites $refSNP_path/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
--known-sites $refSNP_path/1000G_phase1.indels.hg19.sites.vcf \
--known-sites $refSNP_path/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-L $interval_file
#--binary-tag-name ContextCovariate \
#--binary-tag-name CycleCovariate\
# -nct 80 \

#  STEP 4: Apply base quality score recalibration (takes ~3-5 min)
gatk ApplyBQSR \
-I $data_path/${sample}.bam \
-bqsr $output_path/${sample}_recal.table \
--reference $BWA_fasta_ref_file \
-L $interval_file \
-O $output_path/${sample}_recal.bam 
java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
    I=$output_path/${sample}_recal.bam  \
    O=$output_path/${sample}_recal_rg.bam \
    SORT_ORDER=coordinate \
    RGID=psych \
    RGLB=sample \
    RGPU=sample \
    RGPL=illumina \
    RGSM=${sample} 
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_recal_rg.bam



# STEP 5: snp calling for single bam file sample
# 1. Variant  (1-2 min)

gatk HaplotypeCaller  \
--reference $BWA_fasta_ref_file \
-I $output_path/${sample}_recal_rg.bam \
-O $output_path/${sample}.g.vcf.gz \
-L $interval_file \
-ERC GVCF  \
   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
 #
# --annotation-group,-G:String  One or more groups of annotations to apply to variant calls  This argument may be
#                               specified 0 or more times. Default value: null. Possible Values:
#                               {AlleleSpecificAnnotation, AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation,
#                               StandardHCAnnotation, StandardMutectAnnotation}


# 2. Data aggregation step ( for hichip)
# combine gvcfs
 # gatk CombineGVCFs \
 #   -R reference.fasta \
 #   --variant sample1.g.vcf.gz \
 #   --variant sample2.g.vcf.gz \
 #   -O cohort.g.vcf.gz
### -G StandardAnnotation -G AS_StandardAnnotation

# 3. Joint genotyping
# Take the outputs from step 2 (or step 1 if dealing with fewer samples) and run GenotypeGVCFs on all of them together to create the raw SNP and indel VCFs that are usually emitted by the callers.
gatk GenotypeGVCFs \
   -R $BWA_fasta_ref_file \
   -V $output_path/${sample}.g.vcf.gz \
   -O $output_path/${sample}.vcf.gz\
   -G StandardAnnotation -G AS_StandardAnnotation
 

# 4. Variant recalibration
 gatk VariantRecalibrator \
   -R $BWA_fasta_ref_file \
   -V $output_path/${sample}.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${refSNP_path}/hapmap_3.3.hg19.sites.vcf \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ${refSNP_path}/1000G_omni2.5.hg19.sites.vcf \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${refSNP_path}/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${refSNP_path}/dbsnp_138.hg19.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O $output_path/${sample}.recal \
   --tranches-file $output_path/${sample}.tranches \
   --rscript-file $output_path/${sample}.plots.R
 


 gatk ApplyVQSR \
   -R $BWA_fasta_ref_file \
   -V $output_path/${sample}.vcf \
   -O $output_path/${sample}.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file $output_path/${sample}.tranches \
   --recal-file $output_path/${sample}.recal \
   -mode SNP


# -ERC GVCF
# for







