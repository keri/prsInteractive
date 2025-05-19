#!/bin/bash

#
#SBATCH --job-name=prs_score_development
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_score_development.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_score_development.err
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





#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"



module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"


pheno=$1


##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
source ../config.sh  # because you're in prsInteractive/hpc

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
export DATA_TYPE="epi"
python run_model_batches_submit.sh














