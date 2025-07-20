#!/bin/bash

pheno='celiacDisease'


source ../env.config

CONFIG_FILE="$RESULTS_PATH/$pheno/pheno.config"

source $CONFIG_FILE

export HOLDOUT_PATH=$HOLDOUT_PATH
export ENV_FILE=$ENV_FILE
export HLA_FILE=$HLA_FILE
export GENE_ENV_FILE=$GENE_ENV_FILE
export DATA_PATH=$DATA_PATH
export RESULTS_PATH=$RESULTS_PATH

# create hla (and environmental data files?)
#python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"
#
#
##create gene_env datasets for training, validation, and test
#python "$SCRIPTS_DIR/clean_create_environment_data.py"

## Add entries
{
    echo "GENE_ENV_TRAINING=$PHENO_PATH/geneEnvironmentTraining.csv"
    echo "GENE_ENV_TEST=$PHENO_PATH/geneEnvironmentTest.csv"
    echo "GENE_ENV_HOLDOUT=$PHENO_PATH/geneEnvironmentHoldout.csv"
    
} >> "$CONFIG_FILE"

source $CONFIG_FILE

#bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno