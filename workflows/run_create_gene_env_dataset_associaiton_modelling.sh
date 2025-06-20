#!/bin/bash

#
# run_env_gene_creation_association_modelling_submit.sh
#
# This script creates gene_env datasets for training, test, and holdout data



pheno=$1

#PHENO='type2Diabetes_test'
#DATA_TYPE='epi'

echo "PHENO is set to : $pheno"


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

export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export HOLDOUT_PATH=$HOLDOUT_PATH
export PHENO_PATH=$PHENO_PATH
export ENV_FILE=$ENV_FILE
export HLA_FILE=$HLA_FILE
export GENE_ENV_FILE=$GENE_ENV_FILE
#create gene_env datasets for training, validation, and test
python "$SCRIPTS_DIR/clean_create_environment_data.py"


echo "python script for model training started!"
deactivate ukb_env

#run R script for association analysis


