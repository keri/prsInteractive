#!/bin/bash

platform=${5:-"local"}  # default to local
pheno=$1
icd10=$2
phenoStr=$3
n=$4 


#setup environment variables
bash ../envSetUp.sh $pheno $icd10 "${phenoStr}" $n

source ../env.config  # because you're in prsInteractive/workflows

## Set base directories

DATA_DIR="$PROJECT_ROOT/testData"


# Set environment variable
export DATA_PATH="$DATA_DIR"
export PHENO_PATH="$PHENO_DIR"
export RESULTS_PATH="$RESULTS_DIR"
export PHENO="type2Diabetes_test"
export PHENO_STR="type 2 diabetes"
export ICD="E11"
export N_CORES=18


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: $PHENO_STR"



#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_PATH}" ]; then
    echo "Folder '${PHENO_PATH}' does not exist. Creating it..."
    mkdir "${PHENO_PATH}" 

    
else
    echo "Folder '${PHENO_PATH}' already exists."	
fi



###################  SPECIFIC TO TEST #####################

#create genotype data in .bed format
bash "$SCRIPTS_DIR/test/create_simulation_genotype_data.sh"


# create phenotype data and train test split IDs
python3 "$SCRIPTS_DIR/test/create_simulated_participant_covariate_data.py"


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"


# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
# test script to handle special SNPs in chr 6
bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh"

bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"

#only done for test as the IID needs to be numeric and simulation creates IID as string
python "$SCRIPTS_DIR/test/change_test_data_IID_to_number.py"

# Generate conversion file from existing .raw file
# Quick one-liner to fix existing files
for file in "$PHENO_PATH/trainingCombined.raw" "$PHENO_PATH/testCombined.raw" "$PHENO_PATH/holdoutCombined.raw"; do
    awk 'BEGIN{OFS=" "} NR==1{print; next} {gsub(/per/, "", $1); gsub(/per/, "", $2); $1=int($1); $2=int($2); print}' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done


bash run_model_batches_test.sh $PHENO "main"

PHENO_CONFIG="$PHENO_PATH/pheno.config"
source $PHENO_CONFIG
#
EPI_PATH="$PHENO_PATH/epiFiles"
export EPI_PATH
#
python "${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"
#
NEW_EPI_FILE="$PHENO_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary"
#
echo "Epi file is now set to ... $NEW_EPI_FILE"
#
#check to see if EPI_PATH exists
if grep -q "^EPI_FILE=" "$PHENO_CONFIG"; then
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS sed syntax (requires '' for in-place)
        sed -i '' "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$PHENO_CONFIG"
    else
        # Linux sed syntax
        sed -i "s|^EPI_FILE=.*|EPI_FILE=${NEW_EPI_FILE}|" "$PHENO_CONFIG"
    fi
else
    echo "EPI_FILE=${NEW_EPI_PATH}" >> "$PHENO_CONFIG"
    echo "export EPI_FILE" >> "$PHENO_CONFIG"	
fi
#
#
##run the epi batch models on the hpc
export DATA_TYPE="epi"
bash run_model_batches_test.sh

#add 4 lines of epi interactions with significance not present after feature ranking step
python "$SCRIPTS_DIR/test/add_epi_importantFeatures_for_test.py"


#run the GxGxE analysis
export ENV_TYPE='cardioMetabolic'
python "$SCRIPTS_DIR/gene_environment_feature_discovery.py"

#run 
#bash "$SCRIPTS_DIR/run_plink_LD.sh"












