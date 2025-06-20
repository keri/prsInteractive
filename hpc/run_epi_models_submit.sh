#!/bin/bash

#
#SBATCH --job-name=run_epi_models
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_run_epi_models.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_run_epi_models.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=50
#SBATCH --mem=850G
#SBATCH --time=19:00:00
#


pheno=$1


module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"


sbatch run_model_batches_submit.sh $pheno 'epi'










