#!/bin/bash

#
#SBATCH --job-name=prs_score_development
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_score_development.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_score_development.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=150G
#SBATCH --time=3:00:00
#

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/scores/featureWeights.csv
#   $PHENO_PATH/scores/importantFeaturesPostShap.csv
#   $RESULTS_PATH/covar.txt
#   $PHENO_PATH/scores/importantFeaturesPostShap.csv




#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"



module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"

module load plink/1.90

pheno=$1


##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
source ../config.sh  # because you're in prsInteractive/hpc

PHENO_DIR="$RESULTS_DIR/$pheno"

export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"

echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"



###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


if [ ! -d "${PHENO_DIR}" ]; then
    echo "Folder '${PHENO_DIR}' does not exist. The data cleaning and creation step must be complete with run_data_cleaning_workflow..."
    echo "You must know the ICD10 code to filter, the phenotype substring to look for and phenotype string to create results folder (user defined).."
    exit 1
    
else
    echo "Folder '${PHENO_DIR}' already exists."	
fi

#create the GxGxE training, test, and holdout datasets for modeling
#the datasets are created using parameters trained with training set and used to transfrom test and holdout sets
#these are: feature means for E features and imputation before GxE combined, and scaler after GxE combined.
python "${SCRIPTS_DIR}/create_gene_env_data_for_model_training.py"

#r script for running the final models with output being: $PHENO_PATH/scores/finalFeatureWeights.csv
#R finalModelTrainingWithCovariate.R












