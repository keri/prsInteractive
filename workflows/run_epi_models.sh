#!/bin/bash

pheno=$1

# FILES THAT MUST BE PRESENT:
#   $PHENO_PATH/pheno_config.sh
#   $PHENO_PATH/epiFiles/trainingCombinedEpi.epi.cc.summary
#   $PHENO_PATH/trainingCombined.raw
#   $PHENO_PATH/testCombined.raw

##############  SET UP ENV VARIABLES FOR JOB #################

# Source config
source ../config.sh  # because you're in prsInteractive/hpc

PHENO_DIR="$RESULTS_DIR/$pheno"

source "$PHENO_DIR/pheno_config.sh"


#check to see that epi analysis is done and results are present
if [ ! -d "${EPI_PATH}" ]; then
    echo "Folder '${EPI_PATH}' does not exist. Epistatic analysis needs be completed using the multiprocessing_fast_epistasis bash script..."
    exit 1
fi

export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"

echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"

EPI_PATH="$PHENO_PATH/epiFiles"
export EPI_PATH
python "${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"

NEW_EPI_PATH="$PHENO_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary"


echo "Epi path is set to ... NEW_EPI_PATH"
export EPI_PATH=$NEW_EPI_PATH

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
bash run_model_batches.sh














