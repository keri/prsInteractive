#!/bin/bash

pheno=$1

# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/pheno.config
#   $PHENO_PATH/epiFiles/trainingCombinedEpi.epi.cc.summary


##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
source ../env.config # because you're in prsInteractive/hpc

PHENO_DIR="$RESULTS_DIR/$pheno"

source "$PHENO_DIR/pheno.config"

#PHENO_PATH is sourced in pheno.config
EPI_PATH="$PHENO_PATH/epiFiles"
#check to see that epi analysis is done and results are present
if [ ! -d "${EPI_PATH}" ]; then
    echo "Folder '${EPI_PATH}' does not exist. Epistatic analysis needs be completed using the multiprocessing_fast_epistasis bash script..."
    exit 1
fi

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














