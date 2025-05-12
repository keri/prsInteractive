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

#!/bin/bash
#Keri Multerer October 2025
#iterate over batches
#
#SBATCH --job-name=model_batches
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keri@multerer.com
#

source ../config.sh

# Function to display usage information
usage() {
    echo "Usage: $0 [pheno] [data_type] "
    echo "  input_file  : File containing data to be processed in batches and epi or main for data_type"
    echo ""
    echo "Example: $0 type2Diabetes epi"
    exit 1
}

# Check for minimum required arguments
if [ $# -lt 2 ]; then
    usage
fi

PHENO="$1"
DATA_TYPE="$2"


export DATA_TYPE=$DATA_TYPE
export PHENO=$PHENO

if [ "$data_type" == "epi" ]; then
    INPUT_FILE=$EPI_PATH
else
    INPUT_FILE="${PHENO_PATH}/merged_allChromosomes_columns.snplist"
fi


# Validate input file exists 
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi


if [ ! -d "${PHENO_DIR}" ]; then
    echo "${PHENO_DIR} does not exist. You need to go back and create the data by running run_data_cleaning_workflow_submit.sh"
    exit 1
        
else
    echo "Running batch models for ${DATA_TYPE} data... "	
fi

# Calculate total number of lines in the file
TOTAL_LINES=$(wc -l < "$INPUT_FILE")
echo "Total lines in file: $TOTAL_LINES"

# Calculate total number of batches (batch size is 3000)
BATCH_SIZE=3000
TOTAL_BATCHES=$(( (TOTAL_LINES + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "Total number of 3K batches: $TOTAL_BATCHES"

# Loop through batches and submit jobs in groups of 5
BATCHES_PER_JOB=5
TOTAL_JOBS=$(( (TOTAL_BATCHES + BATCHES_PER_JOB - 1) / BATCHES_PER_JOB ))
echo "Grouping into $TOTAL_JOBS jobs (5 batches per job)"

for JOB_ID in $(seq 1 $TOTAL_JOBS); do
    # Calculate batch range for this job
    JOB_START_BATCH=$(( (JOB_ID - 1) * BATCHES_PER_JOB + 1 ))
    JOB_END_BATCH=$(( JOB_ID * BATCHES_PER_JOB ))
    
    # Ensure end batch doesn't exceed total batches
    if [ $JOB_END_BATCH -gt $TOTAL_BATCHES ]; then
        JOB_END_BATCH=$TOTAL_BATCHES
    fi
    
    echo "Submitting job $JOB_ID to process batches $JOB_START_BATCH to $JOB_END_BATCH"
    
    export START=$JOB_START_BATCH
    export END=$JOB_END_BATCH
    # Submit the SLURM job
    sbatch sklearnSectionModelsScoreTrain_submit.sh
#   --export=START_BATCH=$JOB_START_BATCH,END_BATCH=$JOB_END_BATCH,DATA_TYPE=$DATA_TYPE,PHENO=$PHENO \

done

echo "Batch job submission complete!"