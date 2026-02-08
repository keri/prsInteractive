#!/bin/bash

#
#SBATCH --job-name=gene_env_discovery
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_gene_env_discovery.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=80
#SBATCH --mem=500G
#SBATCH --time=85:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#

module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

pheno=$1
INPUT_FILE=$2
CHUNK_START=$3
CHUNK_STOP=$4
EPI_COMBO=${5:-"prod"}
threshold=${6:-1.99}
env_type=${7:-"cardioMetabolic"}



echo "[DEBUG] All required files found. Starting Python script..."

export PHENO_DATA
export TRAINING_PATH
export TEST_PATH
export ENV_TYPE=$env_type
export PHENO=$pheno
export RESULTS_PATH
export WITHDRAWAL_PATH
export CHUNK_START=$CHUNK_START
export CHUNK_STOP=$CHUNK_STOP
export INPUT_FILE=$INPUT_FILE
export THRESHOLD=$threshold
export EPI_COMBO=$EPI_COMBO


    
# Run the Python script
python "${SCRIPTS_DIR}/gene_environment_feature_discovery.py"

exit_code=$?
echo "[DEBUG] Python script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script failed with exit code $exit_code"
    exit $exit_code
fi

echo "[DEBUG] Script started successfully"

    











