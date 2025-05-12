#!/bin/bash

#run the sklearnSectionModelScoreTrain.py with test data

# Source config
source ../testData/config_test.sh  # because you're in prsInteractive/workflows

data_type=$1

pheno='type2Diabetes'

PHENO_DIR="$RESULTS_DIR/$pheno"

export PHENO_PATH="$PHENO_DIR"

echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"

###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


python "${SCRIPTS_PATH}/sklearnSectionModelsScoreTrain.py"