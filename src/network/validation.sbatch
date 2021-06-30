#!/bin/bash
# validation.sbatch
#
#SBATCH -J validation # A single job name for the array
#SBATCH -p normal # Partition
#SBATCH -n 1 # one core
#SBATCH --mem=16000
#SBATCH --gres gpu:1
#SBATCH --time=24:00:00
#SBATCH -o validation.out # Standard output
#SBATCH -e validation.err # Standard error
#SBATCH --mail-type=FAIL,END # notifications for job done & fail
#SBATCH --mail-user=mguo123@stanford.edu


# make & move into new directory, and run!
ml tensorflow/1.1.0-cp36

srun python validation.py .
