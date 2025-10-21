#!/bin/bash

#
#SBATCH --job-name=glmFinalModel
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_shap_values_post_training.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_shap_values_post_training.err
#SBATCH --partition=longrun
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=14:00:00
#


#load required models
module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

pheno=$1

# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi

#check that a results folder for phenotype exists
if [ ! -d "${RESULTS_PATH}/$pheno" ]; then
    echo "Folder '${RESULTS_PATH}/$pheno' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${RESULTS_PATH}/$pheno/pheno.config"
fi


python "${SCRIPTS_DIR}/calculate_shap_values_post_training.py"