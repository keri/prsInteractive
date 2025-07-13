#!/bin/bash

#
# run_model_batches_submit.sh
#
# This script calculates the number of 3K-sized batches needed for processing
# and submits jobs for all batches.
#
# Usage: ./calculate_batches.sh [input_file]
#
# Example: ./calculate_batches.sh my_data.txt

#SBATCH --job-name=model_batches
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keri@multerer.com
#


module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

PHENO=$1
DATA_TYPE=$2
#PHENO='type2Diabetes'
#DATA_TYPE='epi'

echo "PHENO is set to : $PHENO"
echo "DATA_TYPE is set to : $DATA_TYPE"


# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi

#check that a results folder for phenotype exists
if [ ! -d "${RESULTS_PATH}/$pheno" ]; then
    echo "Folder '${RESULTS_PATH}/$pheno' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${RESULTS_PATH}/$PHENO/pheno.config"
fi



echo "Running batch models for $DATA_TYPE data... "	

PHENO_PATH="$RESULTS_PATH/$PHENO"

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
    echo "epi file is set to : $EPI_FILE"
    #check to see if epi file has been filtered for redundant pairs
    if [[ "$EPI_FILE" == *"filtered"* ]]; then
        echo "✓ epi pairs have been filtered for redundancy and file updated"
    else
        #filter redundant epi pairs
        export EPI_FILE=$EPI_FILE
        echo "reading EPI_FILE : ${EPI_FILE}"
        echo "[PYTHON] running python script ${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"
        python "$SCRIPTS_DIR/filter_redundant_epi_pairs.py"
        #wait 10 mins for script to finish
        echo "Python script filtering redundant epi features finished with exit code: $?"
        NEW_EPI_FILE="$EPI_FILE.filtered"
        # Verify filtered file exists
        if [ ! -f "$NEW_EPI_FILE" ]; then
            echo "❌ Expected filtered file not found: $NEW_EPI_FILE"
            echo "Continuing with original EPI file"
        else
            # Update config file
            CONFIG_FILE="${PHENO_PATH}/pheno.config"
            echo "Updating config file: $CONFIG_FILE"
            echo "New EPI_FILE value: $NEW_EPI_FILE"
            
            # Create backup
            cp "$CONFIG_FILE" "${CONFIG_FILE}.backup"
            
            if [[ "$(uname)" == "Darwin" ]]; then
                # macOS version
                sed -i '' "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$CONFIG_FILE"
            else
                # Linux version  
                sed -i "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$CONFIG_FILE"
            fi
            
            # Verify the change
            if grep -q "^EPI_FILE=${NEW_EPI_FILE}$" "$CONFIG_FILE"; then
                echo "✓ Successfully updated EPI_FILE in config"
                rm "${CONFIG_FILE}.backup"  # Remove backup if successful
                INPUT_FILE=$NEW_EPI_FILE
                EPI_FILE=$NEW_EPI_FILE
            else
                echo "❌ Failed to update config file"
                mv "${CONFIG_FILE}.backup" "$CONFIG_FILE"  # Restore backup
                
            fi

        
        fi
        
        
    fi
    
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
BATCH_SIZE=3000
TOTAL_BATCHES=$(( (TOTAL_LINES + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "Total number of 3K batches: $TOTAL_BATCHES"

# Loop through batches and submit jobs in groups of 5
BATCHES_PER_JOB=5
TOTAL_JOBS=$(( (TOTAL_BATCHES + BATCHES_PER_JOB - 1) / BATCHES_PER_JOB ))
echo "Grouping into $TOTAL_JOBS jobs (5 batches per job)"

#or JOB_ID in $(seq 1 $TOTAL_JOBS); do
for JOB_ID in $(seq 1 2); do
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
    export EPI_FILE 
    
    
    # Submit the SLURM job
    #   --export=START_BATCH=$JOB_START_BATCH,END_BATCH=$JOB_END_BATCH,DATA_TYPE=$DATA_TYPE,PHENO=$PHENO \
    sbatch "${HPC_DIR}/sklearnSectionModelsScoreTrain_submit.sh"
    
done

echo "python script for model training started!"