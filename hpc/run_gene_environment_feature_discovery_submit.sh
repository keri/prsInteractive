#!/bin/bash

#
#SBATCH --job-name=prs_score_development
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=00:10:00
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
env_type=$2

##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
source ../env.config  # because you're in prsInteractive/hpc


#export PHENO_PATH="$PHENO_DIR"
export ENV_TYPE=$env_type


echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"

PHENO_CONFIG="$RESULTS_DIR/$pheno/pheno.config.sh"

if [ ! -f "${PHENO_CONFIG}"]; then
    echo "Folder ${PHENO_CONFIG} does not exist. The envSetUp.sh step must be done and data cleaning and creation step completed with run_data_cleaning_workflow..."
    exit 1
    
else
    source "${PHENO_CONFIG}"	
fi


#create the GxGxE training, test, and holdout datasets for modeling
#the datasets are created using parameters trained with training set and used to transfrom test and holdout sets
#these are: feature means for E features and imputation before GxE combined, and scaler after GxE combined.
export PHENO_PATH=$PHENO_PATH
export PHENO=$pheno


echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_DIR"
echo "[HPC WORKFLOW] TRAINING FILE is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST FILE is set: $TEST_PATH"
echo "[HPC WORKFLOW] HOLDOUT FILE is set: $HOLDOUT_PATH"

python "${SCRIPTS_DIR}/gene_environment_feature_discovery.py"













