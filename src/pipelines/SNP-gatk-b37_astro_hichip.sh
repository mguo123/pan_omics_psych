#!/bin/bash
#SNP-gatk-b37-hichip_Astrocytes.sbatch
#SBATCH --job-name=gatk-hichip-{cell_type}
#SBATCH --output=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/hichip_output_snp/Astrocytes_gatk.out
#SBATCH --error=/oak/stanford/groups/khavari/users/mguo123/psych_project/data/hichip_output_snp/Astrocytes_gatk.err
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

BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/b37/human_g1k_v37.fasta"
# BWA_fasta_ref_file="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg38/Homo_sapiens_assembly38.fasta"
# BWA_fasta_ref_dir="/oak/stanford/groups/khavari/users/mguo123/reference_genome/hg19_UCSC"
BWA_fasta_ref_dir="/oak/stanford/groups/khavari/users/mguo123/reference_genome/b37"
picard_path="/oak/stanford/groups/rbaltman/mguo123/picard/build/libs"

refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/b37"

# refSNP_path="/oak/stanford/groups/khavari/users/mguo123/refsnps/hg19"
# # make a merged hichip bed file of all the tissue anchor peaks (from FitHiChIP so a bit stringent but that's okay) see OneNote for more details
interval_file="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/hicpro_longrange/all_neuro_pks_b37.bed"
# consider using only atac peak regions???
# interval_file="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/atac_bed_files/all_neuro_pks_b37.bed"




# SET SAMPLE/celltype INFO
# # Celltypes: Astrocytes 
cell_type="Astrocytes"
sample_list=(Astro_B1 Astro_B2)
data_path="/oak/stanford/groups/khavari/users/lkhd/project/3D/HiChIP/novogene/hicpro/outputs/astro_output/bowtie_results/bwt2"


### STEP 1A: setting output path
output_path="/oak/stanford/groups/khavari/users/mguo123/psych_project/data/hichip_output_snp"
output_path=${output_path}/$cell_type
echo $output_path
mkdir -p $output_path

cd $output_path



for sample in ${sample_list[@]}; do 
   ## FOR EVERYTHING EXCEPT ASTROCYTES
# sample=${sample_list[${idx}]}
echo $sample

### STEP 2 merge the bam files per sample of a tissue
bam_files=`ls ${data_path}/${sample}/*bwt2pairs.bam`

merge_cmd="java -jar ${picard_path}/picard.jar MergeSamFiles O=${output_path}/${sample}_merged.bam"
for f in $bam_files; do 
echo $f
merge_cmd="${merge_cmd} I=${f}"
done
echo $merge_cmd
eval $merge_cmd
echo 'done merging bam files of sample: ${sample}...'
 
# # STEP 2 BAM PREPROCESSING (TAKES A WHILE) make bam files and index (deduped)
# ## alignment (takes 90-hours)
# #http://bio-bwa.sourceforge.net/bwa.shtml
# bwa mem -R "@RG\tID:psych\tSM:sample\tPL:illumina\tLB:sample\tPU:sample" ${BWA_fasta_ref_dir}/hg19.fa ${fastqList_1} ${fastqList_2} > $output_path/${sample}.sam
# # bwa mem -R "@RG\tID:psych\tSM:sample\tPL:illumina\tLB:sample\tPU:sample" ${BWA_fasta_ref_dir}/hg19.fa /oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/H9iN-day0-A-GCTACGCT_S3_L004_R1_001.fastq.gz /oak/stanford/groups/khavari/users/lkhd/project/Wernig/ATAC/rawdata/H9iN-day0-A-GCTACGCT_S3_L004_R2_001.fastq.gz > $output_path/${sample}.sam

# add read group, dedup (already done), sort index (~15 min)
# java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
#     I=$bam_file  \
#     O=$output_path/${sample}.bam \
#     SORT_ORDER=coordinate \
#     RGID=psych \
#     RGLB=sample \
#     RGPU=sample \
#     RGPL=illumina \
#     RGSM=${sample} 
# java -jar $picard_path/picard.jar SortSam INPUT=$output_path/${sample}.bam OUTPUT=$output_path/${sample}_sorted.bam SORT_ORDER=coordinate

java -jar $picard_path/picard.jar SortSam INPUT=${output_path}/${sample}_merged.bam OUTPUT=$output_path/${sample}_sorted.bam SORT_ORDER=coordinate
java -jar $picard_path/picard.jar MarkDuplicates \
      I=$output_path/${sample}_sorted.bam \
      O=$output_path/${sample}_sorted_mkdup.bam \
      M=$output_path/${sample}_marked_dup_metrics.txt
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_sorted_mkdup.bam

java -jar $picard_path/picard.jar AddOrReplaceReadGroups \
    I=$output_path/${sample}_sorted_mkdup.bam   \
    O=$output_path/${sample}_sorted_mkdup_rg.bam  \
    SORT_ORDER=coordinate \
    RGID=psych \
    RGLB=sample \
    RGPU=sample \
    RGPL=illumina \
    RGSM=${sample}
java -jar $picard_path/picard.jar BuildBamIndex INPUT=$output_path/${sample}_sorted_mkdup_rg.bam



## DEBUG CHECK Readgroups
# samtools view -H $output_path/${sample}_sorted_mkdup.bam  | grep '@RG'




# STEP 3: base recalibration takes 5-10 min
#First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
gatk BaseRecalibrator \
-R $BWA_fasta_ref_file \
-I $output_path/${sample}_sorted_mkdup_rg.bam \
-O $output_path/${sample}_recal.table \
--known-sites $refSNP_path/dbsnp_138.b37.vcf \
--known-sites $refSNP_path/1000G_phase1.snps.high_confidence.b37.vcf \
--known-sites $refSNP_path/1000G_phase1.indels.b37.vcf \
--known-sites $refSNP_path/Mills_and_1000G_gold_standard.indels.b37.vcf \
-L $interval_file


#  STEP 4: Apply base quality score recalibration (takes ~3-5 min)
gatk ApplyBQSR \
-I $output_path/${sample}_sorted_mkdup_rg.bam \
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

#  ### FOR THE NGN2 just copy (or if only 1 rep)
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
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${refSNP_path}/hapmap_3.3.b37.vcf \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ${refSNP_path}/1000G_phase1.indels.b37.vcf \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${refSNP_path}/1000G_phase1.snps.high_confidence.b37.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${refSNP_path}/dbsnp_138.b37.vcf \
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


