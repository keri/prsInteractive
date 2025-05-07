#!/bin/bash



# Set base directories
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$WORKFLOW_DIR")"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
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
python3 "$SCRIPTS_DIR/create_pheno_train_test_split.py"


# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
# test script to handle special SNPs in chr 6
bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh"

bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"








