#!/bin/bash


PHENO=$1

#PHENO='type2Diabetes'

EPI_COMBO=${2:-"sum"}
#EPI_COMBO="prod"

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
    "$SCRIPTS_DIR/prsAUCDelongStats.R"
    "${SCRIPTS_DIR}/calculate_top_features_in_cohort.py"
    "${SCRIPTS_DIR}/calculate_prs_stats.py"
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


source "${CONFIG_FILE}"

#check to see if filtered file is present
if [ ! -f "$PHENO_DATA/scores/featureScoresReducedFinalModel.filtered.csv" ]; then
    echo "[G-E non-additive FILTERING] is not done"
    echo "run filter_non_additive_gene_env_features.py"
    #filter non-additive gene-env features
    python "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py" \
    --config_file "$PHENO_DATA/pheno.config" \
    --feature_scores_file "$PHENO_DATA/scores/featureScoresReducedFinalModel.csv"
else
    echo "✓ GxGxE features have been filtered for non-additive post modelling"
fi

##check to see if LD has been done previously before association
#if [ ! -f "$PHENO_PATH/finalModel.ld" ];then
#   bash "$SCRIPTS_DIR/run_plink_LD.sh" $PHENO
#fi
#
##necessary exports are done at top of script
##Run the Python script
#
python "${SCRIPTS_DIR}/calculate_prs_post_modelling.py"
#
#exit_code=$?
#echo "[DEBUG] Python script calculating prs exited with code: $exit_code"
#
#if [ $exit_code -ne 0 ]; then
#   echo "ERROR: Python script calculating prs failed with exit code $exit_code"
#   exit $exit_code
#fi
#
##export of PHENO_DATA done

python "${SCRIPTS_DIR}/combine_prs.py"
--pheno_data $PHENO_DATA

exit_code=$?
echo "[DEBUG] Python combine_prs.py script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python combine_prs.py script failed with exit code $exit_code"
    exit $exit_code
fi

#calculate statstics: McNemar, ttest, pearson correlation, precision/recall improvement over G and exclusive cases 
#calculate peformance of trained models and PRS calculations for main v other
#need PHENO_DATA exported or added as an --pheno_data argument
python "${SCRIPTS_DIR}/calculate_prs_stats.py" \
--pheno_data $PHENO_DATA

RScript "$SCRIPTS_DIR/prsAUCDelongStats.R" \
--pheno_data $PHENO_DATA

########  CALCULATE TOP FEATURES IN STATISTICALLY DISTINCT COHORTS ##########

python "${SCRIPTS_DIR}/calculate_top_features_in_cohort.py" \
--pheno_data $PHENO_DATA \
--feature_scores_file_filtered "$PHENO_DATA/scores/featureScoresReducedFinalModel.filtered.csv" \
--threshold 1.99


echo "[DEBUG] Script completed successfully"












