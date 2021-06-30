#!/bin/bash
#
#SBATCH --job-name=bigwig
#SBATCH --output=/oak/stanford/groups/khavari/users/lkhd/project/Vanessa/RNA/bigwig/sample.out
#SBATCH --error=/oak/stanford/groups/khavari/users/lkhd/project/Vanessa/RNA/bigwig/sample.err
#SBATCH --time=10:00:00
#SBATCH -p khavari,owners,normal
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=lkhd@stanford.edu

ml biology
ml bedtools

for sam in IRF6-NT1-A IRF6-NT1-B IRF6-sg1-A IRF6-sg1-B IRF6-sg2-A IRF6-sg2-B IRF6-sg8-A IRF6-sg8-B Low-1A Low-1B Low-3A Low-3B Norm-1A Norm-1B Norm-3A Norm-3B;
do
dir=/oak/stanford/groups/khavari/users/lkhd/project/Vanessa/RNA/STAR_genome/
out=/oak/stanford/groups/khavari/users/lkhd/project/Vanessa/RNA/bigwig/
mkdir ${out}${sam}
cd ${out}${sam}
genomeCoverageBed -bg -ibam ${dir}${sam}/${sam}_Aligned.sortedByCoord.out.bam -g /oak/stanford/groups/khavari/users/namyoung/anno/ucsc/hg19.chrom.sizes > ${sam}.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ${sam}.bedGraph > ${sam}_sorted_rna.bedGraph
/oak/stanford/groups/khavari/users/lkhd/software/bedGraphToBigWig ${sam}_sorted_rna.bedGraph /oak/stanford/groups/khavari/users/namyoung/anno/STAR/R19/indice-100/chrNameLength.txt ${sam}.rna.bw
done
