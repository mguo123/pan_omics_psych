#!/bin/bash
# protein_protein_rand_cmd.sbatch
#
#SBATCH -J PPI_rand # A single job name for the array
#SBATCH -p PPI_rand # Partition
#SBATCH -n 1 # one core
#SBATCH --mem=16000
#SBATCH --gres gpu:1
#SBATCH --time=24:00:00
#SBATCH -o protein_rand_col%A_%a.out # Standard output
#SBATCH -e protein_rand_col%A_%a.err # Standard error
#SBATCH --mail-type=FAIL,END # notifications for job done & fail
#SBATCH --mail-user=mguo123@stanford.edu


# the inputted number is the number of random seed nodes to run
ml python/2.7.5
printf -v SLURM_ARRAY_TASK_ID "%02d" $SLURM_ARRAY_TASK_ID ;
echo ${SLURM_ARRAY_TASK_ID}
# mkdir -p outputs_scratch
python protein_protein_cmd.py DATA/new_edges_t0.5_binary DATA/seed_nodes.csv 0 1 ${SLURM_ARRAY_TASK_ID}

# to run: sbatch --array=0-71 protein_protein_cmd.sbatch