#!/bin/bash

#
#SBATCH --job-name=calculate_prs
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_prs.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_prs.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=oconnoas@staff.vuw.ac.nz
#



module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

pheno=$1


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
export TEST_PATH=$TEST_PATH
export HOLDOUT_PATH=$HOLDOUT_PATH
export PHENO=$pheno
export RESULTS_PATH=$RESULTS_PATH
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH
export DATA_PATH=$DATA_PATH
export HLA_FILE=$HLA_FILE
export COVAR_FILE=$COVAR_FILE
export TEST_ENV_GEN_FILE=$GENE_ENV_TEST
export HOLDOUT_ENV_GEN_FILE=$GENE_ENV_HOLDOUT
export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_PATH: $PHENO_PATH"
echo "SCRIPTS_DIR: $SCRIPTS_DIR"
echo "TEST_PATH: $TEST_PATH"
echo "HOLDOUT_PATH: $HOLDOUT_PATH"
echo "RESULTS_PATH: $RESULTS_PATH"
echo "TEST_ENV_GEN_FILE: $TEST_ENV_GEN_FILE"
echo "HOLDOUT_ENV_GEN_FILE: $HOLDOUT_ENV_GEN_FILE"
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
    "${SCRIPTS_DIR}/calculate_prs_for_filtered_main_epi.py"
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
python "${SCRIPTS_DIR}/calculate_prs_for_filtered_main_epi.py"

exit_code=$?
echo "[DEBUG] Python script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script failed with exit code $exit_code"
    exit $exit_code
fi

python "${SCRIPTS_DIR}/combine_prs.py" 

exit_code=$?
echo "[DEBUG] Python combine_prs.py script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python combine_prs.py script failed with exit code $exit_code"
    exit $exit_code
fi

echo "[DEBUG] Script completed successfully"












