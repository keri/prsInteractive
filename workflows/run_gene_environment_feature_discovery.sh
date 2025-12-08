#!/bin/bash



pheno=$1
env_type=$2
threshold=${3:-2}
EPI_COMBO=${4:-"sum"}

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


#check that a results folder for phenotype exists
if [ ! -d "${PHENO_DATA}" ]; then
    echo "Folder '${PHENO_DATA}' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${PHENO_DATA}/pheno.config"
fi


echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"


# Set phenotype-specific paths


echo "[WORKFLOW] ENV_TYPE is set to: $env_type"
export PHENO_PATH="$RESULTS_PATH/$pheno"
export PHENO_DATA=$PHENO_DATA
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export ENV_TYPE=$env_type
export PHENO=$pheno
export RESULTS_PATH=$RESULTS_PATH


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_DATA: $PHENO_DATA"
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
    "$PHENO_DATA/scores/importantFeaturesPostShap.csv"
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

export INPUT_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.csv"
export THRESHOLD=$threshold
eport EPI_COMBO
# Run the Python script
python "${SCRIPTS_DIR}/gene_environment_feature_discovery.py"

exit_code=$?
echo "[DEBUG] Python script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script failed with exit code $exit_code"
    exit $exit_code
fi

echo "[DEBUG] Script completed successfully"











