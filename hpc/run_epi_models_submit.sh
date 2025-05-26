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
#   $PHENO_PATH/pheno_config.sh
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
if [ ! -f "../config.sh" ]; then
    echo "ERROR: ../config.sh not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../config.sh 2>/dev/null || echo '../config.sh')"
    exit 1
fi

source ../config.sh


PHENO_DIR="$RESULTS_DIR/$pheno"

PHENO_CONFIG="$PHENO_DIR/pheno_config.sh"

source "$PHENO_CONFIG"

export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"

echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"


###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


if [ ! -d "${PHENO_DIR}" ]; then
    echo "Folder '${PHENO_DIR}' does not exist. The data cleaning and creation step must be complete with run_data_cleaning_workflow..."
    echo "You must know the ICD10 code to filter, the phenotype substring to look for and phenotype string to create results folder (user defined).."
    exit 1
    
else
    echo "Folder '${PHENO_DIR}' already exists."	
fi

#filter the epi pairs in .summary file to remove redundanct pairs that are identical but reversed

EPI_PATH="$PHENO_PATH/epiFiles"
export EPI_PATH

python "${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"

TIMEOUT=600

conda deactivate

TIMEOUT=15

NEW_EPI_PATH="$PHENO_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary"
export EPI_PATH=$NEW_EPI_PATH
echo "Epi path is now set to ... $NEW_EPI_PATH"

#check to see if EPI_PATH exists
if grep -q "^EPI_PATH=" "$PHENO_CONFIG"; then
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS sed syntax (requires '' for in-place)
        sed -i '' "s|^EPI_PATH=.*|EPI_PATH=${NEW_EPI_PATH}|" "$PHENO_CONFIG"
    else
        # Linux sed syntax
        sed -i "s|^EPI_PATH=.*|EPI_PATH=${NEW_EPI_PATH}|" "$PHENO_CONFIG"
    fi
else
    echo "EPI_PATH=${NEW_EPI_PATH}" >> "$PHENO_CONFIG"
    echo "export EPI_PATH" >> "$PHENO_CONFIG"	
fi


#run the epi batch models on the hpc

sbatch run_model_batches_submit.sh $PHENO 'epi'













