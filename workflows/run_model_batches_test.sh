#!/bin/bash

# run_model_batches_submit.sh
#
# This script calculates the number of 3K-sized batches needed for processing
# and submits jobs for all batches.



# Function to display usage information
#usage() {
#   echo "Usage: $0 [pheno] [data_type] "
#   echo "  input_file  : File containing data to be processed in batches and epi or main for data_type"
#   echo ""
#   echo "Example: $0 type2Diabetes epi"
#   #exit 1
#}
#
## Check for minimum required arguments
#if [ $# -lt 2 ]; then
#   usage
#fi

PHENO="$1"
echo "PHENO is set to : $PHENO"
DATA_TYPE="$2"
echo "DATA_TYPE is set to : $DATA_TYPE"

PHENO_CONFIG="$PHENO_PATH/pheno.config"

# Validate pheno_config file exists and source if it does
if [ ! -f $PHENO_CONFIG ]; then
    echo "$PHENO_CONFIG does not exist which means pheno.config has not been produced ... "
    echo "You need to back and run run_data_cleaning_workflow_submit.sh "
    exit 1
else 
    #exports PHENO_PATH, EPI_PATH, TRAINING_PATH, TEST_PATH
    source $PHENO_CONFIG
fi

#check that pheno path exists
if [ ! -d "$PHENO_PATH" ]; then
    echo "'$PHENO_PATH' does not exist. You need to go back and create the data by running run_data_cleaning_workflow_submit.sh"
    exit 1
    
else
    echo "Running batch models for '${DATA_TYPE}' data... "	
fi

# Validate model folders exist in pheno folder
if [ ! -d "$PHENO_PATH/scores" ]; then
    echo "creating '$PHENO_PATH/scores' folder ... "
    mkdir "$PHENO_PATH/scores"
fi

if [ ! -d "$PHENO_PATH/models" ]; then
    echo "creating '$PHENO_PATH/models' folder ... "
    mkdir "$PHENO_PATH/models"
fi

if [ ! -d "$PHENO_PATH/figures" ]; then
    echo "creating '$PHENO_PATH/figures' folder ... "
    mkdir "$PHENO_PATH/figures"
fi


if [ "$DATA_TYPE" == "epi" ]; then
    INPUT_FILE=$EPI_FILE
    echo "in the run batch models: EPI_PATH is set to: $EPI_FILE"
else
    INPUT_FILE="$PHENO_PATH/merged_allChromosomes.snplist"
fi


# Validate input file exists 
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

TOTAL_LINES=$(wc -l < "$INPUT_FILE")

echo "Total lines: $TOTAL_LINES"

# Calculate total number of batches (batch size is 3000)
BATCH_SIZE=$TOTAL_LINES
TOTAL_BATCHES=$(( (TOTAL_LINES + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "Total number of batches: $TOTAL_BATCHES"

# Loop through batches and submit jobs in groups of 5
BATCHES_PER_JOB=1
TOTAL_JOBS=$(( (TOTAL_BATCHES + BATCHES_PER_JOB - 1) / BATCHES_PER_JOB ))
echo "Grouping into $TOTAL_JOBS jobs (5 batches per job)"
    
#for JOB_ID in $(seq 1 $TOTAL_BATCHES); do
for JOB_ID in $(seq 1 1); do
    echo "job id : $JOB_ID"
    # Calculate batch range for this job
    JOB_START_BATCH=$(( (JOB_ID - 1) * BATCHES_PER_JOB + 1 ))
    JOB_END_BATCH=$(( JOB_ID * BATCHES_PER_JOB ))
    
    echo "job ending in batch : $JOB_END_BATCH"
    # Ensure end batch doesn't exceed total batches
    if [ $JOB_END_BATCH -gt $TOTAL_BATCHES ]; then
        JOB_END_BATCH=$TOTAL_BATCHES
    fi
    
    echo "Submitting job $JOB_ID to process batches $JOB_START_BATCH to $JOB_END_BATCH"
    
    export START=$JOB_START_BATCH
    export END=$JOB_END_BATCH
    export DATA_TYPE
    export PHENO
    export PHENO_PATH
    export TRAINING_PATH
    export TEST_PATH


    # Submit the SLURM job
    #   --export=START_BATCH=$JOB_START_BATCH,END_BATCH=$JOB_END_BATCH,DATA_TYPE=$DATA_TYPE,PHENO=$PHENO \
    python "${SCRIPTS_DIR}/sklearnSectionModelsScoreTrain.py"
    
done
    
echo "Batch job submission complete!"