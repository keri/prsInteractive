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
env_type=$2

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



# Check if RESULTS_PATH was set from env.config
if [ -z "$RESULTS_PATH" ]; then
    echo "ERROR: RESULTS_PATH not set from env.config"
    exit 1
fi

echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"


# Set phenotype-specific paths


echo "[WORKFLOW] ENV_TYPE is set to: $env_type"
export PHENO_PATH="$RESULTS_PATH/$pheno"
export ENV_TYPE=$env_type
export PHENO=$pheno


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_PATH: $PHENO_PATH"
echo "ENV_TYPE: $ENV_TYPE"
echo "====================================="



INPUT_FILE="$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv"
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