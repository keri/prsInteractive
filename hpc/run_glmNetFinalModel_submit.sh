#!/bin/bash

#
#SBATCH --job-name=glmFinalModel
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_glmFinalModel.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_glmFinalModel.err
#SBATCH --partition=longrun
#SBATCH --cpus-per-task=80
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#


#load required models
module load GCC/11.2.0  
module load OpenMPI/4.1.1
module load R/4.2.0

module load plink/1.90

pheno=$1
# Source config
source ../env.config # because you're in prsInteractive/hpc


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

export RESULTS_PATH=$RESULTS_PATH
export PHENO_PATH=$PHENO_PATH
export DATA_PATH=$DATA_PATH
export HLA_FILE="$RESULTS_PATH/participant_hla.csv"
export COVAR_FILE="$RESULTS_PATH/covar.csv"
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export TRAINING_ENV_GEN_FILE=$GENE_ENV_TRAINING
export TEST_ENV_GEN_FILE=$GENE_ENV_TEST

bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno

# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--results_path "$RESULTS_PATH" \
--pheno_path "$PHENO_PATH" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" 


# Add entries
{
    echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_SCORES=$PHENO_PATH/scores/modelScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scores/predictProbsReducedFinalModel.csv"
} >> "${RESULTS_PATH}/$pheno/pheno.config"


bash "run_prs_calculations_submit.sh" $pheno