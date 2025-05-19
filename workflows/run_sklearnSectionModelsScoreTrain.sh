#!/bin/bash

#run the sklearnSectionModelScoreTrain.py with test data

# Source config
echo "[HPC WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "starting iteration ...: $START"
echo "ending iteration ...: $END"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_DIR"
echo "[HPC WORKFLOW] TRAINING_PATH is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST_PATH is set: $TEST_PATH"


###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


python "${SCRIPTS_DIR}/sklearnSectionModelsScoreTrain.py"