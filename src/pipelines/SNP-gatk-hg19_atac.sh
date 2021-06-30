#!/bin/bash
#SNP-gatk-hg19-atac.sbatch
#SBATCH --job-name=gatk-atac-{cell_type}
#SBATCH --output=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp/${cell_type}_gatk.out
#SBATCH --error=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp/${cell_type}.err
#SBATCH --time=48:00:00
#SBATCH -p rbaltman,khavari,normal,owners
#SBATCH --nodes=1
#SBATCH --mem=100000
#SBATCH -c 16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mguo1238@stanford.edu
# source ~/.bashrc
module load biology
ml gatk
ml bwa
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
# BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg38/Homo_sapiens_assembly38.fasta"
BWA_fasta_ref_dir="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19_UCSC"
 # /oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19

refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/hg19"
# # make a merged atac bed file of all the tissue atac peaks
# /oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files]$ cat H9iN-day0.narrowPeak.bed H9iN-day28.narrowPeak.bed H9-Ngn2.narrowPeak.bed SLC_baseline.narrowPeak.bed SLC-Ngn2.narrowPeak.bed SL-Ngn2.narrowPeak.bed H9_baseline.narrowPeak.bed H9iN-day10.narrowPeak.bed H9iN-day4.narrowPeak.bed SL_baseline.narrowPeak.bed  SLC.narrowPeak.bed SL.narrowPeak.bed | bedtools sort -i stdin| bedtools merge -i stdin > all_neuro_pks.bed
interval_file="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files/all_neuro_pks.bed"
picard_path="/oak/stanford/groups/rbaltman/mguo123/picard/build/libs"

## SET SAMPLE/celltype INFO
# Celltypes: Astrocytes H9_D0 H9_D10 H9_D28 H9_D2 H9_D4 SLC_D0 SLC_D2 SL_D0 SL_D2 
# note the H9_D0 done in sdev 
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(H9iN-day0-A H9iN-day0-B)
cell_type="H9_D0"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(H9iN-day10-A H9iN-day10-B)
cell_type="H9_D10"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(H9iN-day28-A H9iN-day28-B)
cell_type="H9_D28"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(H9iN-day4-A H9iN-day4-B)
cell_type="H9_D4"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(SL-A SL-B)
cell_type="SL"
data_path="/oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/"
sample_list=(SLC-A SLC-B)
cell_type="SLC"


data_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/cheen_data_ngn2/atac_data/starting_fastq/"
sample_list=(H9-Ngn2)
cell_type="H9_D2"
data_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/cheen_data_ngn2/atac_data/starting_fastq/"
sample_list=(SL-Ngn2)
cell_type="SL_D2"
data_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/cheen_data_ngn2/atac_data/starting_fastq/"
sample_list=(SLC-Ngn2)
cell_type="SLC_D2"

data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/ATAC/sfgf/hiseq_051519/"
sample_list=(Astrocytes-B1 Astrocytes-B2)
cell_type="Astrocytes"






### STEP 1A: setting output path
output_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_output_snp"
output_path=${output_path}/$cell_type
echo $output_path
mkdir -p $output_path



for sample in ${sample_list[@]}; do echo $sample; 

# STEP 1B
# list fastq files
cd $data_path

fastqList_1=$(find "$(pwd)" -name "${sample}*_R1_001.fastq.gz" )
fastqList_2=$(find "$(pwd)" -name "${sample}*_R2_001.fastq.gz" )

# ### FOR ASTROCYTES ONLY (and thus all epithelial cells for cancer project)
fastqList_1=$(find "$(pwd)" -name "${sample}*[0-9]_merged_L001_R1_001.fastq.gz" )
fastqList_2=$(find "$(pwd)" -name "${sample}*[0-9]_merged_L001_R2_001.fastq.gz" )

echo $fastqList_1
echo $fastqList_2
cd $output_path


# # STEP 2 BAM PREPROCESSING (TAKES A WHILE) make bam files and index (deduped)
## alignment (takes 90-hours)
#http://bio-bwa.sourceforge.net/bwa.shtml
bwa mem -R "@RG\tID:psych\tSM:sample\tPL:illumina\tLB:sample\tPU:sample" ${BWA_fasta_ref_dir}/hg19.fa ${fastqList_1} ${fastqList_2} > $output_path/${sample}.sam
# bwa mem -R "@RG\tID:psych\tSM:sample\tPL:illumina\tLB:sample\tPU:sample" ${BWA_fasta_ref_dir}/hg19.fa /oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/H9iN-day0-A-GCTACGCT_S3_L004_R1_001.fastq.gz /oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/H9iN-day0-A-GCTACGCT_S3_L004_R2_001.fastq.gz > $output_path/${sample}.sam

# add read group, dedup (already done), sort index (~15 min)

java -jar $picard_path/picard.jar SortSam INPUT=$output_path/${sample}.sam OUTPUT=$output_path/${sample}_sorted.bam SORT_ORDER=coordinate
java -jar $picard_path/picard.jar MarkDuplicates \
      I=$output_path/${sample}_sorted.bam \
      O=$output_path/${sample}_sorted_mkdup.bam \
      M=$output_path/${sample}_marked_dup_metrics.txt
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_sorted_mkdup.bam



## DEBUG CHECK Readgroups
# samtools view -H $output_path/${sample}_sorted_mkdup.bam  | grep '@RG'




# STEP 3: base recalibration takes 5-10 min
#First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
gatk BaseRecalibrator \
-R $BWA_fasta_ref_file \
-I $output_path/${sample}_sorted_mkdup.bam \
-O $output_path/${sample}_recal.table \
--known-sites $refSNP_path/dbsnp_138.hg19.vcf \
--known-sites $refSNP_path/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
--known-sites $refSNP_path/1000G_phase1.indels.hg19.sites.vcf \
--known-sites $refSNP_path/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-L $interval_file


#  STEP 4: Apply base quality score recalibration (takes ~3-5 min)
gatk ApplyBQSR \
-I $output_path/${sample}_sorted_mkdup.bam \
-bqsr $output_path/${sample}_recal.table \
--reference $BWA_fasta_ref_file \
-L $interval_file \
-O $output_path/${sample}_recal.bam 
# java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
#     I=$output_path/${sample}_recal.bam  \
#     O=$output_path/${sample}_recal_rg.bam \
#     SORT_ORDER=coordinate \
#     RGID=psych \
#     RGLB=sample \
#     RGPU=sample \
#     RGPL=illumina \
#     RGSM=${sample} 
# java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_recal_rg.bam



# STEP 5: snp calling for single bam file sample
# allele specific version : https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants
# 1. Variant  (20 min) get gvcf
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531812?id=4017 - GVCF vs VCF

gatk HaplotypeCaller  \
--reference $BWA_fasta_ref_file \
-I $output_path/${sample}_recal.bam \
-O $output_path/${sample}.g.vcf.gz \
-L $interval_file \
-ERC GVCF  \
   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
 #
# --annotation-group,-G:String  One or more groups of annotations to apply to variant calls  This argument may be
#                               specified 0 or more times. Default value: null. Possible Values:
#                               {AlleleSpecificAnnotation, AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation,
#                               StandardHCAnnotation, StandardMutectAnnotation}


done




# STEP6. Data aggregation step ( for hichip)
# combine gvcfs
gatk CombineGVCFs \
--reference $BWA_fasta_ref_file \
--variant $output_path/${sample_list[0]}.g.vcf.gz \
--variant $output_path/${sample_list[1]}.g.vcf.gz \
-O $output_path/${cell_type}.g.vcf.gz
 -G StandardAnnotation -G AS_StandardAnnotation

#  ### FOR THE NGN2 just copy
# mv $output_path/${sample_list[0]}.g.vcf.gz $output_path/${cell_type}.g.vcf.gz

# STEP7. Joint genotyping
# Take the outputs from step 2 (or step 1 if dealing with fewer samples) and run GenotypeGVCFs on all of them together to create the raw SNP and indel VCFs that are usually emitted by the callers.
gatk GenotypeGVCFs \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.g.vcf.gz \
   -O $output_path/${cell_type}.vcf.gz\
   -G StandardAnnotation -G AS_StandardAnnotation
 

# STEP 8 . Variant recalibration
# vcf spec https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
# https://gatk.broadinstitute.org/hc/en-us/articles/360037594511-VariantRecalibrator#top
 gatk VariantRecalibrator \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.vcf.gz \
   -AS \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${refSNP_path}/hapmap_3.3.hg19.sites.vcf \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ${refSNP_path}/1000G_omni2.5.hg19.sites.vcf \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${refSNP_path}/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${refSNP_path}/dbsnp_138.hg19.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -mode SNP \
   -O $output_path/${cell_type}.AS.recal \
   --tranches-file $output_path/${cell_type}.AS.tranches \
   --rscript-file $output_path/${cell_type}.plots.AS.R
 
 gatk ApplyVQSR \
   -R $BWA_fasta_ref_file \
   -V $output_path/${cell_type}.vcf.gz \
   -O $output_path/${cell_type}_postvqsr.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file $output_path/${cell_type}.AS.tranches \
   --recal-file $output_path/${cell_type}.AS.recal \
   -mode SNP








