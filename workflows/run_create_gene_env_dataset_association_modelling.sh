#!/bin/bash

#pheno=$1
#env_type=$2
#threshold=${3:-2}
threshold=1.99
pheno='celiacDisease'
env_type='cardioMetabolic'
EPI_COMBO='prod'
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
export GENE_ENV_FILE="$PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShap.csv"
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH
export THRESHOLD=$threshold
export EPI_COMBO=$EPI_COMBO

if [ -f "$GENE_ENV_FILE" ]; then
    echo "[DEBUG] Filtered features file found"
    
    {
        echo "GENE_ENV_FILE=$GENE_ENV_FILE"
    } >> "$CONFIG_FILE"
else
    echo "ERROR: Filtered features file not found at $GENE_ENV_FILE"
    exit 1
fi

#filter features using threshold
python "$SCRIPTS_DIR/find_important_gene_environment_featuress.py" --gene_env_file $GENE_ENV_FILE --threshold 2

#check to see if GxGxE features have been filtered by Shap value
if [ -f "$PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShapFilteredZscore.csv" ]; then
    echo "[DEBUG] Filtered features file found"
    
    {
        echo "GENE_ENV_FILE_FILTERED_ZSCORE=$PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShapFilteredZscore.csv"
    } >> "$CONFIG_FILE"
    
else
    echo "ERROR: Filtered features file not found at $PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShapFilteredZscore.csv"
    echol "Run filtered GxGxE using a threshold"
    exit 1
fi

export GENE_ENV_FILE_FILTERED_ZSCORE=$GENE_ENV_FILE_FILTERED_ZSCORE
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
    echo "CLINICAL_ENV_TRAINING=$PHENO_DATA/clinicalEnvironmentTraining.csv"
    echo "CLINICAL_ENV_TEST=$PHENO_DATA/clinicalEnvironmentTest.csv"
    echo "CLINICAL_ENV_HOLDOUT=$PHENO_DATA/clinicalEnvironmentHoldout.csv"

} >> "$CONFIG_FILE"



echo "python script for model training started!"
#deactivate ukb_env

#run R script for association analysis
