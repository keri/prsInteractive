#!/bin/bash

#pheno=$1
#env_type=$2
#threshold=${3:-2}
threshold=2
pheno='celiacDisease'
env_type='cardioMetabolic'
EPI_COMBO='sum'
#EPI_COMBO=${4:-"sum"}

# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi


if [ "${EPI_COMBO}" == "sum" ]; then
    COMBO_FOLDER='summedEpi'
else
    COMBO_FOLDER='productEpi'
fi

PHENO_DATA="${RESULTS_PATH}/$pheno/${COMBO_FOLDER}"
CONFIG_FILE="${PHENO_DATA}/pheno.config"

#check that a results folder for phenotype exists
if [ ! -d "${PHENO_DATA}" ]; then
    echo "Folder '${PHENO_DATA}' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${CONFIG_FILE}"
fi


echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "PHENO is set to : $pheno"


export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export HOLDOUT_PATH=$HOLDOUT_PATH
export PHENO_DATA=$PHENO_DATA
export ENV_FILE=$ENV_FILE
export HLA_FILE=$HLA_FILE
export GENE_ENV_FILE=$GENE_ENV_FILE
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH
export THRESHOLD=$threshold

#create gene_env datasets for training, validation, and test
python "$SCRIPTS_DIR/clean_create_environment_data.py"

## Add entries to CONFIG_FILE
{
    echo "GENE_ENV_TRAINING=$PHENO_DATA/geneEnvironmentTraining.csv"
    echo "GENE_ENV_TEST=$PHENO_DATA/geneEnvironmentTest.csv"
    echo "GENE_ENV_HOLDOUT=$PHENO_DATA/geneEnvironmentHoldout.csv"
    echo "MAIN_ENV_TRAINING=$PHENO_DATA/mainEnvironmentTraining.csv"
    echo "MAIN_ENV_TEST=$PHENO_DATA/mainEnvironmentTest.csv"
    echo "MAIN_ENV_HOLDOUT=$PHENO_DATA/mainEnvironmentHoldout.csv"
    echo "ALL_ENV_TRAINING=$PHENO_DATA/allEnvironmentTraining.csv"
    echo "ALL_ENV_TEST=$PHENO_DATA/allEnvironmentTest.csv"
    echo "ALL_ENV_HOLDOUT=$PHENO_DATA/allEnvironmentHoldout.csv"
    
} >> "$CONFIG_FILE"


#update config file with filtered env data
## Add entries
if [ -f "$PHENO_DATA/cardioMetabolicimportantFeaturesPostShapFiltered.csv" ]; then
    echo "[DEBUG] Filtered features file found"
    
    {
        echo "FILTERED_GENE_ENV_FEATURES=$PHENO_DATA/cardioMetabolicimportantFeaturesPostShapFiltered.csv"
    } >> "$CONFIG_FILE"
else
    echo "ERROR: Filtered features file not found at $PHENO_DATA/cardioMetabolicimportantFeaturesPostShapFiltered.csv"
    exit 1
fi


echo "python script for model training started!"
#deactivate ukb_env

#run R script for association analysis
