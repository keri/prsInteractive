#!/bin/bash

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   raw variant calls with bed format : ukb{project_number}_c1_b0_v2
#   raw hla data in csv format  : hla_participant.csv',index_col='Participant ID') and file with headers : ukb_hla_v2.txt
#   participant data in csv format with 

#pheno=$1
#icd10=$2
#phenoStr=$3
#n=$4

#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"
#n=40

#source ../config.sh #because you're in workflows directory
set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Source environment if it exists
if [ -f ".env" ]; then
    source .env
    echo "Loaded environment from .env"
else
    echo "No .env file found, running setup..."
    ./setup-environment.sh
    source .env
fi

# Set default cromwell jar location if not specified
CROMWELL_JAR="${CROMWELL_JAR:-$HOME/cromwell/cromwell-90.jar}"

# Check if cromwell jar exists
if [ ! -f "$CROMWELL_JAR" ]; then
    echo "Error: Cromwell JAR not found at $CROMWELL_JAR"
    echo "Please set CROMWELL_JAR environment variable or place cromwell jar at the default location"
    exit 1
fi

# Export the project root for Docker mounting
export PRS_INTERACTIVE_HOME="$SCRIPT_DIR"

echo "Running prsInteractive pipeline..."
echo "Project root: $PRS_INTERACTIVE_HOME"
echo "Cromwell JAR: $CROMWELL_JAR"

# Run the workflow
java -jar "$CROMWELL_JAR" run \
workflows/prsInteractivePipeline.DataProcessing.wdl \
-i workflows/pipelineInputs.json \
-o workflows/options.json \
--conf workflows/cromwell.conf

echo "Pipeline execution completed!"


# Set base directories
PHENO_DIR="$RESULTS_DIR/$pheno"

export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"
export PHENO_STR="$phenoStr"
export ICD="$icd10"


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: $PHENO_STR"



#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_DIR}" ]; then
    echo "Folder '${PHENO_DIR}' does not exist. Creating it..."
    mkdir "${PHENO_DIR}" 

    
else
    echo "Folder '${PHENO_DIR}' already exists."	
fi


#create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"

# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
bash "$SCRIPTS_DIR/plink_clean_variant_calls.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh"

bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"

export DATA_TYPE="main"
bash run_model_batches.sh









