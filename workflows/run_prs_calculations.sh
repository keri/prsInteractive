#!/bin/bash


#PHENO=$1

PHENO='celiacDisease'

EPI_COMBO=${2:-"sum"}

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

PHENO_DATA="${RESULTS_PATH}/$PHENO/${COMBO_FOLDER}"
CONFIG_FILE="${PHENO_DATA}/pheno.config"

#check that a results folder for phenotype exists
if [ ! -d "${PHENO_DATA}" ]; then
    echo "Folder '${PHENO_DATA}' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $PHENO env variables."
    #source pheno specific environment variables
    source "${CONFIG_FILE}"
fi

echo "[WORKFLOW]PHENO_DATA is to : $PHENO_DATA"
echo "[WORKFLOW]TRAINING_PATH is to : $TRAINING_PATH"
echo "[WORKFLOW]TEST_PATH is to : $TEST_PATH"
echo "[WORKFLOW]PHENO is to : $PHENO"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"


# Set phenotype-specific paths


echo "[WORKFLOW] ENV_TYPE is set to: $env_type"
export PHENO_DATA
export PHENO_PATH
export TEST_PATH
export HOLDOUT_PATH
export PHENO
export RESULTS_PATH
export WITHDRAWAL_PATH
export DATA_PATH
export HLA_FILE
export COVAR_FILE
export GENE_ENV_TEST
export GENE_ENV_HOLDOUT
export GENE_MAIN_TEST
export GENE_MAIN_HOLDOUT
export FEATURE_SCORES_FILE
export ENV_FILE
export EPI_COMBO


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_DATA: $PHENO_DATA"
echo "SCRIPTS_DIR: $SCRIPTS_DIR"
echo "TEST_PATH: $TEST_PATH"
echo "HOLDOUT_PATH: $HOLDOUT_PATH"
echo "RESULTS_PATH: $RESULTS_PATH"
echo "TEST_ENV_GEN_FILE: $GENE_ENV_TEST"
echo "HOLDOUT_ENV_GEN_FILE: $GENE_ENV_HOLDOUT"
echo "TEST_ENV_MAIN_FILE: $MAIN_ENV_TEST"
echo "HOLDOUT_ENV_MAIN_FILE: $MAIN_ENV_HOLDOUT"
echo "FEATURE_SCORES_FILE : $FEATURE_SCORES_FILE"
echo "====================================="

# Check if required files exist before running Python
echo "[DEBUG] Checking required files..."

required_files=(
    "$FEATURE_SCORES_FILE"
    "$GENE_ENV_HOLDOUT"
    "$GENE_ENV_TEST"
    "$WITHDRAWAL_PATH"
    "$COVAR_FILE"
    "$HLA_FILE"
    "$HOLDOUT_PATH"
    "$TEST_PATH"
    "${SCRIPTS_DIR}/calculate_prs_post_modelling.py"
    "$SCRIPTS_DIR/run_plink_LD.sh"
    "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file does not exist: $file"
        exit 1
    else
        echo "[DEBUG] Found: $file"
    fi
done

echo "[DEBUG] All required files found. Starting Python scripts..."

#check if file exists post GxGxE filtering
for file in "$PHENO_DATA/scores/featureScoresReducedFinalModel.filtered.csv"; do
    if [[ -f "$file" ]]; then
        echo "✓ File exists: $file"
        
    else
        echo "✗ File missing: $file"
        export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
        export CONFIG_FILE=$CONFIG_FILE
        python "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py"
    fi
done

source "${CONFIG_FILE}"

#
##### ENSURE CONFIG FILE IS UPDATED ##########
if [[ "$FEATURE_SCORES_FILE" == *"filtered"* ]]; then
    echo "✓ GxGxE features have been filtered for non-additive post modelling"
else
    #wait 10 mins for script to finish
    echo "Python script filtering non-additive GxGxE features finished with exit code: $?"
    echo "Continuing with original FEATURE SCORES file"
fi

#check to see if LD has been done previously before association
if [ ! -f "$PHENO_DATA/finalModel.ld" ];then
    export PHENO_DATA=$PHENO_DATA
    export PHENO_PATH=$PHENO_PATH
    bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
fi

#necessary exports are done at top of script
#Run the Python script

python "${SCRIPTS_DIR}/calculate_prs_post_modelling.py"

exit_code=$?
echo "[DEBUG] Python script calculating prs exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script calculating prs failed with exit code $exit_code"
    exit $exit_code
fi

#export of PHENO_DATA done
python "${SCRIPTS_DIR}/combine_prs.py"

exit_code=$?
echo "[DEBUG] Python combine_prs.py script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python combine_prs.py script failed with exit code $exit_code"
    exit $exit_code
fi

#calculate statstics: McNemar, ttest, pearson correlation, precision/recall improvement over G and exclusive cases 
python "${SCRIPTS_DIR}/calculate_prs_stats.py"

exit_code=$?
echo "[DEBUG] Python calculate_prs_stats.py script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python calculate_prs_stats.py script failed with exit code $exit_code"
    exit $exit_code
fi

#calculate AUC Delong statistics
Rscript "$SCRIPTS_DIR/prsAUCDelongStats.R" \
--pheno_data "$PHENO_DATA" 

exit_code=$?
echo "[DEBUG] Rscript prsAUCDelongStats.R script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Rscript prsAUCDelongStats.R script failed with exit code $exit_code"
    exit $exit_code
fi

########  CALCULATE TOP FEATURES IN STATISTICALLY DISTINCT COHORTS ##########
export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
python "${SCRIPTS_DIR}/calculate_top_features_in_cohort.py" 


echo "[DEBUG] Script completed successfully"












