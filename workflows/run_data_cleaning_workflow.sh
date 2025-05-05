#!/bin/bash

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   raw variant calls with bed format : ukb{project_number}_c1_b0_v2
#   raw hla data in csv format  : hla_participant.csv',index_col='Participant ID') and file with headers : ukb_hla_v2.txt
#   participant data in csv format with 

pheno=$1
icd10=$2
phenoStr=$3


# Set base directories
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$WORKFLOW_DIR")"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR = "$PROJECT_ROOT/results"
PHENO_DIR="$RESULTS_DIR/$pheno/"


# Set environment variable
export DATA_PATH="$DATA_DIR"
export RESULTS_PATH="$RESULTS_DIR"
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
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: $PHENO"



#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_PATH}" ]; then
    echo "Folder '${PHENO_PATH}' does not exist. Creating it..."
    mkdir "${PHENO_PATH}" 

    
else
    echo "Folder '${PHENO_PATH}' already exists."	
fi


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"

# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environmental_hla_covar_data.py"

# Run the variant call cleaning
bash "$SCRIPTS_DIR/plink_clean_variant_calls.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes_submit.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/plink_convert_merged_to_A_submit.sh"

#run fast epistasis on merged bed files
bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis_submit.sh"

# create filtered snp list
#bash "$SCRIPTS_DIR/create_filtered_snp_list_2.sh"







