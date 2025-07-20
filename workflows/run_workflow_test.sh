#!/bin/bash

platform=${5:-"local"}  # default to local
pheno="type2Diabetes_test"
icd10="E11"
phenoStr="type 2 diabetes"
n=18

#setup environment variables
bash ../envSetUp.sh $pheno $icd10 "${phenoStr}" $n

source ../env.config  # because you're in prsInteractive/workflows


#
### Set base directories
mkdir -p "$PRS_INTERACTIVE_HOME/testData/variant_calls"
chmod +x "$PRS_INTERACTIVE_HOME/testData"
DATA_PATH="$PRS_INTERACTIVE_HOME/testData"
PHENO_PATH=$PRS_INTERACTIVE_HOME/results/$pheno
WITHDRAWAL_PATH=/Users/kerimulterer/prsInteractive/testData/withdrawals.csv


#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_PATH}" ]; then
	echo "Folder '${PHENO_PATH}' does not exist. Creating it..."
	mkdir "${PHENO_PATH}" 

else
	echo "Folder '${PHENO_PATH}' already exists."	
fi

PHENO_CONFIG="$PHENO_PATH/pheno.config"
if [ ! -f "${PHENO_CONFIG}" ]; then
	echo "FILE '${PHENO_CONFIG}' does not exist. You need to ensure ../envSetUp.sh has been run with appropriate arguments"
	
else
	source "${PHENO_CONFIG}"
fi

# Set environment variable
export PRS_INTERACTIVE_HOME=$PRS_INTERACTIVE_HOME
export DATA_PATH="$DATA_PATH"
export PHENO_PATH="$PHENO_PATH"
export RESULTS_PATH="$RESULTS_PATH"
export PHENO=$pheno
export PHENO_STR="${phenoStr}"
export ICD=$icd10
export N_CORES=18
export SCRIPTS_DIR=$SCRIPTS_DIR
export EPI_PATH=$EPI_PATH
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH

#
echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: ${PHENO_STR}"


###################  SPECIFIC TO TEST #####################

#create genotype data in .bed format
bash "$SCRIPTS_DIR/test/create_simulation_genotype_data.sh"


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/test/create_simulated_participant_covariate_data.py"


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"


# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"


CONFIG_FILE="$PHENO_PATH/pheno.config"

# Add entries
{
	echo "HLA_FILE=$RESULTS_PATH/participant_hla.csv"
	echo "COVAR_FILE=$RESULTS_PATH/covar.csv" 
	echo "ENV_FILE=$RESULTS_PATH/participant_environment.csv"
} >> "$CONFIG_FILE"

# Simple verification
echo "Config file updated. Contents:"
cat "$CONFIG_FILE"

# Check if actual files exist
for file in "$RESULTS_PATH/participant_hla.csv" "$RESULTS_PATH/covar.csv" "$RESULTS_PATH/participant_environment.csv"; do
	if [[ -f "$file" ]]; then
		echo "✓ File exists: $file"
	else
		echo "✗ File missing: $file"
	fi
done

# Run the variant call cleaning
# test script to handle special SNPs in chr 6
bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh"

#source PHENO_CONFIG AGAIN TO GET THE EPI_PATH
bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"


#only done for test as the IID needs to be numeric and simulation creates IID as string
python "$SCRIPTS_DIR/test/change_test_data_IID_to_number.py"


# Quick one-liner to fix existing files
for file in "$PHENO_PATH/trainingCombined.raw" "$PHENO_PATH/testCombined.raw" "$PHENO_PATH/holdoutCombined.raw"; do
	awk 'BEGIN{OFS=" "} NR==1{print; next} {gsub(/per/, "", $1); gsub(/per/, "", $2); $1=int($1); $2=int($2); print}' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done


bash run_model_batches.sh $PHENO "main"


#run the epi batch models on the hpc
bash run_model_batches.sh $PHENO 'epi'


#Validate pheno_config file exists 
if [ ! -f "$PHENO_PATH/pheno.config" ]; then
	echo " '$PHENO_PATH/pheno.config' is not present re run previous steps... "
	
else
	source "$PHENO_PATH/pheno.config"
fi


#run the GxGxE analysis
export ENV_TYPE='cardioMetabolic'
export TEST_PATH=$TEST_PATH
export TRAINING_PATH=$TRAINING_PATH
python "$SCRIPTS_DIR/gene_environment_feature_discovery.py"

#update pheno.config
CONFIG_FILE="$PHENO_PATH/pheno.config"

# Add entries
{
	echo "GENE_ENV_FILE=$PHENO_PATH/scores/cardioMetabolicimportantFeaturesPostShap.csv"
	echo "ENV_TYPE=$ENV_TYPE"
} >> "$CONFIG_FILE"

source $CONFIG_FILE

export HOLDOUT_PATH=$HOLDOUT_PATH
export ENV_FILE=$ENV_FILE
export HLA_FILE=$HLA_FILE
export GENE_ENV_FILE=$GENE_ENV_FILE
#create gene_env datasets for training, validation, and test
python "$SCRIPTS_DIR/clean_create_environment_data.py"


## Add entries
{
	echo "GENE_ENV_TRAINING=$PHENO_PATH/geneEnvironmentTraining.csv"
	echo "GENE_ENV_TEST=$PHENO_PATH/geneEnvironmentTest.csv"
	echo "GENE_ENV_HOLDOUT=$PHENO_PATH/geneEnvironmentHoldout.csv"

} >> "$CONFIG_FILE"

source $CONFIG_FILE

bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno

export PHENO_PATH=$PHENO_PATH
# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--results_path "$RESULTS_PATH" \
--pheno_path "$PHENO_PATH" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" 


# Add entries
{
	echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModel.csv"
	echo "FINAL_MODEL_SCORES=$PHENO_PATH/scoresmodelScoresReducedFinalModel.csv"
	echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scorespredictProbsReducedFinalModel.csv"
} >> "$CONFIG_FILE"

source $CONFIG_FILE

#calculate PRS from model association and after LD filter
export HLA_FILE=$HLA_FILE
export COVAR_FILE=$COVAR_FILE
export GENE_ENV_TEST=$GENE_ENV_TEST
export GENE_ENV_HOLDOUT=$GENE_ENV_HOLDOUT
export TEST_PATH=$TEST_PATH
export HOLDOUT_PATH=$HOLDOUT_PATH
export PHENO_PATH=$PHENO_PATH
export PHENO=$PHENO
export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
export SCORES_PATH="$PHENO_PATH/scores"

python "${SCRIPTS_DIR}/calculate_prs_for_filtered_main_epi.py"


conda deactivate

echo "DONE!!"










    