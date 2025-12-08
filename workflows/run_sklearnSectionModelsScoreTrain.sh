#!/bin/bash

#run the sklearnSectionModelScoreTrain.py with test data

# Source config
echo "[HPC WORKFLOW] PHENO_DATA is set to: $PHENO_DATA"
echo "starting iteration ...: $START"
echo "ending iteration ...: $END"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_DIR"
echo "[HPC WORKFLOW] TRAINING_PATH is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST_PATH is set: $TEST_PATH"
echo "[HPC WORKFLOW] WITHDRAWAL_PATH is set to $WITHDRAWAL_PATH"
echo "[HPC WORKFLOW] DATA_TYPE being modelled is $DATA_TYPE"


###########  export env variables for use in python script ############

export PHENO_DATA
export START
export END
export TEST_PATH
export TRAINING_PATH
export WITHDRAWAL_PATH
export DATA_TYPE

if [ "${DATA_TYPE}" == "main" ];then
	EPI_FILE="none"
else
	echo "[HPC WORKFLOW] EPI_FILE being read is $EPI_FILE"
fi

export EPI_FILE
python "${SCRIPTS_DIR}/sklearnSectionModelsScoreTrain.py"