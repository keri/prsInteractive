#!/bin/bash

#
#SBATCH --job-name=gene_env_discovery
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=60
#SBATCH --mem=900G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=oconnoas@staff.vuw.ac.nz
#

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/scores/importantFeaturesPostShap.csv
#   $RESULTS_PATH/testCombined.raw
#   $RESULTS_PATH/trainingCombined.raw
#   /env.config
#   $PHENO_PATH/pheno.config



module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

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
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export ENV_TYPE=$env_type
export PHENO=$pheno
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













