#!/bin/bash

source ../config.sh  # because you're in prsInteractive/workflows

## Set base directories

DATA_DIR="$PROJECT_ROOT/testData"
RESULTS_DIR="$PROJECT_ROOT/testResults"
PHENO_DIR="$RESULTS_DIR/type2Diabetes"

# Set environment variable
export DATA_PATH="$DATA_DIR"
export PHENO_PATH="$PHENO_DIR"
export RESULTS_PATH="$RESULTS_DIR"
export PHENO="type2Diabetes"
export PHENO_STR="type 2 diabetes"
export ICD="E11"
export N=18


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
#bash "$SCRIPTS_DIR/test/create_simulation_genotype_data.sh"


# create phenotype data and train test split IDs
#python3 "$SCRIPTS_DIR/test/create_simulated_participant_covariate_data.py"


# create phenotype data and train test split IDs
#python "$SCRIPTS_DIR/create_pheno_train_test_split.py"


# create hla (and environmental data files?)
#python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
# test script to handle special SNPs in chr 6
#bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

#merge separate chromosome files into one
#bash "$SCRIPTS_DIR/merge_chromosomes.sh"

#bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"

#only done for test as the IID needs to be numeric and simulation creates IID as string
#python "$SCRIPTS_DIR/test/change_test_data_IID_to_number.py"

#export DATA_TYPE="main"
#bash run_model_batches_test.sh 

PHENO_CONFIG="$PHENO_PATH/pheno_config.sh"
source $PHENO_CONFIG

EPI_PATH="$PHENO_PATH/epiFiles"
export EPI_PATH

python "${SCRIPTS_DIR}/filter_redundant_epi_pairs.py"

NEW_EPI_PATH="$PHENO_PATH/epiFiles/trainingCombinedEpi.filtered.epi.cc.summary"

echo "Epi path is now set to ... $NEW_EPI_PATH"

#check to see if EPI_PATH exists
if grep -q "^EPI_PATH=" "$PHENO_CONFIG"; then
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS sed syntax (requires '' for in-place)
        sed -i '' "s|^EPI_PATH=.*|EPI_PATH=${NEW_EPI_PATH}|" "$PHENO_CONFIG"
    else
        # Linux sed syntax
        sed -i "s|^EPI_PATH=.*|EPI_PATH=${NEW_EPI_PATH}|" "$PHENO_CONFIG"
    fi
else
    echo "EPI_PATH=${NEW_EPI_PATH}" >> "$PHENO_CONFIG"
    echo "export EPI_PATH" >> "$PHENO_CONFIG"	
fi


#run the epi batch models on the hpc
export DATA_TYPE="epi"
bash run_model_batches_test.sh

#add 4 lines of epi interactions with significance not present after feature ranking step
python "$SCRIPTS_DIR/test/add_epi_importantFeatures_for_test.py"


#run the GxGxE analysis
export ENV_TYPE='cardioMetabolic'
python "$SCRIPTS_DIR/gene_environment_feature_discovery.py"

#run 
bash "$SCRIPTS_DIR/run_plink_LD.sh"










