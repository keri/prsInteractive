#!/bin/bash

#
#SBATCH --job-name=gene_env_discovery
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=00:10:00
#SBATCH --mail-type=BEGIN,END,FAIL
#

pheno=$1
env_type=${2:-"cardioMetabolic"}
threshold=${3:-1.99}
EPI_COMBO=${4:-"prod"}

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


echo "[WORKFLOW] ENV_TYPE is set to: $env_type"
export PHENO_PATH="$RESULTS_PATH/$pheno"
export PHENO_DATA=$PHENO_DATA
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export ENV_TYPE=$env_type
export PHENO=$pheno
export RESULTS_PATH=$RESULTS_PATH
export ENV_FILE=$ENV_FILE


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


INPUT_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.csv"
echo "input file is : $INPUT_FILE"

TOTAL_LINES=$(wc -l < "$INPUT_FILE")

echo "Total features: $TOTAL_LINES"

# Calculate total number of batches (batch size is 3000)
BATCH_SIZE=400
TOTAL_BATCHES=$(( (TOTAL_LINES + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "Total number of $BATCH_SIZE feature batches: $TOTAL_BATCHES"


echo "[DEBUG] All required files found. Starting Python script..."


for JOB_ID in $(seq 0 $((TOTAL_BATCHES - 1))); do
    echo "job id : $JOB_ID"
    
    echo "Submitting job $JOB_ID to process "
    
    CHUNK_START=$((JOB_ID * BATCH_SIZE))
    CHUNK_STOP=$((CHUNK_START + BATCH_SIZE))
    
    # Run the Python script
    sbatch run_gene_environment_feature_discovery_submit.sh $pheno $env_type $CHUNK_START $CHUNK_STOP $INPUT_FILE
    
    exit_code=$?
    echo "[DEBUG] Python script exited with code: $exit_code"
    
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Python script failed with exit code $exit_code"
        exit $exit_code
    fi
    
done