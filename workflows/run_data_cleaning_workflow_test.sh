#!/bin/bash

pheno=$1
icd10=$2
phenoStr=$3


# Set base directories
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$WORKFLOW_DIR")"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="$PROJECT_ROOT/testData"
PHENO_DIR="$PROJECT_ROOT/testData/$pheno/"


# Set environment variable
export DATA_PATH="$DATA_DIR"
export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"
export PHENO_STR="$phenoStr"
export ICD="$icd10"


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
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
python "$SCRIPTS_DIR/test/create_simulated_participant_covariate_data.py"

bash "$SCRIPTS_DIR/test/create_simulation_genotype_data.sh"


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"

# Run the variant call cleaning
# test script to handle special SNPs in chr 6
bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes_submit.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/plink_convert_merged_to_A_submit.sh"

#run fast epistasis on merged bed files
bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis_submit.sh"

# create filtered snp list
#bash "$SCRIPTS_DIR/create_filtered_snp_list_2.sh"







