#!/bin/bash

#
#SBATCH --job-name=epi_model_run
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=00:15:00
#

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/pheno.config
#   $PHENO_PATH/epiFiles/trainingCombinedEpi.epi.cc.summary
#   $PHENO_PATH/trainingCombined.raw
#   $PHENO_PATH/testCombined.raw


# Make sure we're in the correct directory when the script runs
# This ensures relative paths work correctly regardless of where sbatch is called from

module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"


pheno=$1


##############  SET UP ENV VARIABLES FOR JOB #################

# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi

#RESULTS_DIR is sourced from env.config
PHENO_DIR="$RESULTS_DIR/$pheno"

source "$PHENO_DIR/pheno.config"

#PHENO_PATH is sourced in pheno.config
EPI_PATH="$PHENO_PATH/epiFiles"
#check to see that epi analysis is done and results are present
if [ ! -d "${EPI_PATH}" ]; then
    echo "Folder '${EPI_PATH}' does not exist. Epistatic analysis needs be completed using the multiprocessing_fast_epistasis bash script..."
    exit 1
fi

#sourced from pheno.config
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"

#check to see if epiPairs have been filtered 
if [[ "$EPI_FILE" == *"filtered"* ]]; then
    echo "✓ epi pairs have been filtered for redundancy and file updated"
else
    #filter redundant epi pairs
    python "${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"
    NEW_EPI_FILE="$PHENO_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary"
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS sed syntax (requires '' for in-place)
        sed -i '' "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$PHENO_DIR/pheno.config"
    else
        # Linux sed syntax
        sed -i "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$PHENO_DIR/pheno.config"
    fi
fi

export EPI_FILE
export PHENO
#run the epi batch models on the hpc
bash run_model_batches.sh $PHENO 'epi'

conda deactivate



TIMEOUT = 600







