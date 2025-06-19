#!/bin/bash

#run_gene_environment_discovery.sh

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/scores/importantFeaturesPostShap.csv
#   $RESULTS_PATH/covar.txt
#   $RESULTS_PATH/testCombined.raw
#   $RESULTS_PATH/trainingCombined.raw
#   $RESULTS_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary
#   /env.config
#   $PHENO_PATH/pheno.config



pheno=$1
env_type=$2

echo "[DEBUG] Starting script with pheno=$pheno, env_type=$env_type"

##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
echo "[DEBUG] Sourcing ../env.config"
source ../env.config



# Check if RESULTS_PATH was set from env.config
if [ -z "$RESULTS_PATH" ]; then
    echo "ERROR: RESULTS_PATH not set from env.config"
    exit 1
fi

echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"

# FIXED: Use RESULTS_PATH consistently
PHENO_CONFIG="$RESULTS_PATH/$pheno/pheno.config"
echo "[DEBUG] Looking for pheno config at: $PHENO_CONFIG"

if [ ! -f "${PHENO_CONFIG}" ]; then
    echo "ERROR: File ${PHENO_CONFIG} does not exist. The envSetUp.sh step must be done and data cleaning and creation step completed with run_data_cleaning_workflow..."
    exit 1
else
    echo "[DEBUG] Sourcing pheno config: $PHENO_CONFIG"
    source "${PHENO_CONFIG}"    
fi

# Set phenotype-specific paths

export ENV_TYPE=$env_type
echo "[WORKFLOW] ENV_TYPE is set to: $ENV_TYPE"
export PHENO_PATH="$RESULTS_PATH/$pheno"
export PHENO=$pheno
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export RESULTS_PATH=$RESULTS_PATH


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_PATH: $PHENO_PATH"
echo "SCRIPTS_DIR: $SCRIPTS_DIR"
echo "TRAINING_PATH: $TRAINING_PATH"
echo "TEST_PATH: $TEST_PATH"
echo "HOLDOUT_PATH: $HOLDOUT_PATH"
echo "RESULTS_PATH: $RESULTS_PATH"
echo "ENV_TYPE: $ENV_TYPE"
echo "====================================="

# Check if required files exist before running Python
echo "[DEBUG] Checking required files..."

required_files=(
    "$TRAINING_PATH"
    "$TEST_PATH"
    "$PHENO_PATH/scores/importantFeaturesPostShap.csv"
    "${SCRIPTS_DIR}/gene_environment_feature_discovery.py"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file does not exist: $file"
        exit 1
    else
        echo "[DEBUG] Found: $file"
    fi
done

echo "[DEBUG] All required files found. Starting Python script..."

# Run the Python script
python "${SCRIPTS_DIR}/gene_environment_feature_discovery.py"

exit_code=$?
echo "[DEBUG] Python script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script failed with exit code $exit_code"
    exit $exit_code
fi

echo "[DEBUG] Script completed successfully"













