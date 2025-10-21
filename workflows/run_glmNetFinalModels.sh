#!/bin/bash


platform=${2:-"local"} 
#pheno=$1
pheno='type2Diabetes'
# Source config
source ../env.config # because you're in prsInteractive/hpc

# Generate a fixed configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

## Export environment variable for Docker mounting
#export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"
#
#echo "[SETUP] Setting up environment for phenotype: $pheno"
#echo "[SETUP] Platform: $platform"
#
## Initialize conda for bash (if not already done)
#if ! command -v conda &> /dev/null; then
#   echo "ERROR: conda not found. Please install conda/miniconda first."
#   exit 1
#fi
#
## Initialize conda in this shell
#eval "$(conda shell.bash hook)"
#
## Check if environment exists, create if it doesn't
#ENV_NAME="ukb_env"
#if ! conda env list | grep -q "^${ENV_NAME}\s"; then
#   echo "[CONDA] Creating conda environment from environment.yml..."
#   conda env create -f "$PROJECT_ROOT/environment.yml"
#else
#   echo "[CONDA] Environment $ENV_NAME already exists"
#   # Optionally update environment
#   echo "[CONDA] Updating environment from environment.yml..."
#   conda env update -f "$PROJECT_ROOT/environment.yml" --prune
#   # Activate the environment
#   echo "[CONDA] Activating environment: $ENV_NAME"
#   conda activate $ENV_NAME
#   
#   # Verify key tools are available
#   echo "[CONDA] Verifying environment setup..."
#   python --version
#   R --version
#   echo "[CONDA] Environment activated successfully"
#fi

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



### CHECK TO SEE IF LD HAS BEEN DONE
#if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
#   #run LD script
#   export PRE_POST_ASSOCIATION='pre'
#   bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
#fi
#
##check to see LD was completed
#if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
#   FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesPostShap.csv"
#else
#   FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv"
#fi
#
## Pass to R script
#Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
#--data_path "$DATA_PATH" \
#--hla_file "$HLA_FILE" \
#--covar_file "$COVAR_FILE" \
#--results_path "$RESULTS_PATH" \
#--pheno_path "$PHENO_PATH" \
#--training_env_gen_file "$GENE_ENV_TRAINING" \
#--test_env_gen_file "$GENE_ENV_TEST" \
#--training_file "$TRAINING_PATH" \
#--test_file "$TEST_PATH" 
#--feature_model_file "$FEATURES_FOR_MODEL_FILE"

# Add entries
#{
#   echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModel.csv"
#   echo "FINAL_MODEL_SCORES=$PHENO_PATH/scores/modelScoresReducedFinalModel.csv"
#   echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scores/predictProbsReducedFinalModel.csv"
#} >> "${RESULTS_PATH}/$pheno/pheno.config"
    
bash run_prs_calculations.sh $pheno

#calculate peformance of trained models and PRS calculations for main v other
#need PHENO_PATH exported or added as an --pheno_path argument
python "${SCRIPTS_DIR}/calculate_prs_stats.py"