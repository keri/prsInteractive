#!/bin/bash

platform=${6:-"local"}  # default to local
pheno="type2Diabetes_test"
icd10="E11"
phenoStr="type 2 diabetes"
n=18
EPI_COMBO=${5:-"sum"}

#setup environment variables
bash ../scripts/test/envSetUpTest.sh $pheno $icd10 "${phenoStr}" $n $EPI_COMBO

source ../env.config  # because you're in prsInteractive/workflows


if [ "${EPI_COMBO}" == "sum" ]; then
	COMBO_FOLDER='summedEpi'
else
	COMBO_FOLDER='productEpi'
fi
PHENO_DATA="${RESULTS_PATH}/$pheno/${COMBO_FOLDER}"

CONFIG_FILE="$PHENO_DATA/pheno.config"
if [ ! -f "${CONFIG_FILE}" ]; then
	echo "FILE '${CONFIG_FILE}' does not exist. You need to ensure ../envSetUpTest.sh has been run with appropriate arguments"
	
else
	source "${CONFIG_FILE}"
fi

#
### Set base directories
#mkdir -p "$PRS_INTERACTIVE_HOME/testData/variant_calls"
#chmod +x "$PRS_INTERACTIVE_HOME/testData"
#DATA_PATH="$PRS_INTERACTIVE_HOME/testData"
#PHENO_PATH=$PRS_INTERACTIVE_HOME/results/$pheno
#WITHDRAWAL_PATH=/Users/kerimulterer/prsInteractive/testData/withdrawals.csv


#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_PATH}" ]; then
	echo "Folder '${PHENO_PATH}' does not exist. Run bash ../scripts/test/envSetUpTest.sh to create it..."
	#mkdir "${PHENO_PATH}" 

else
	echo "Folder '${PHENO_PATH}' already exists."	
fi


# Set environment variable
export PRS_INTERACTIVE_HOME
export DATA_PATH
export PHENO_PATH
export RESULTS_PATH
export PHENO=$pheno
export PHENO_STR="${phenoStr}"
export ICD=$icd10
export N_CORES
export SCRIPTS_DIR
export EPI_PATH
export WITHDRAWAL_PATH
export PHENO_DATA

#
echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: ${PHENO_STR}"


###################  SPECIFIC TO TEST #####################

#
#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh" $pheno

#source PHENO_CONFIG AGAIN TO GET THE EPI_PATH
bash "$SCRIPTS_DIR/multiprocessing_fast_epistasis.sh"

#source new environment variable created for EPI_FILE
source "${CONFIG_FILE}"


bash run_model_batches.sh $PHENO "main"


#run the epi batch models on the hpc
bash run_model_batches.sh $PHENO 'epi'

#before GxGxE feature discovery - ensure file is present from batch modelling


echo "[DEBUG] All required files found. Starting Python script..."

#run the GxGxE analysis
ENV_TYPE='cardioMetabolic'
#all other env variables are exported previously
threshold=1
bash run_gene_environment_feature_discovery.sh $PHENO $ENV_TYPE $threshold

# Add entries
{
	echo "GENE_ENV_FILE=$PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShap.csv"
	echo "ENV_TYPE=$ENV_TYPE"
} >> "$CONFIG_FILE"


#update pheno.config
#CONFIG_FILE="$PHENO_DATA/pheno.config"

source $CONFIG_FILE

export HOLDOUT_PATH
export ENV_FILE
export HLA_FILE
export GENE_ENV_FILE
#create gene_env datasets for training, validation, and test
bash run_create_gene_env_dataset_association_modelling.sh $PHENO $ENV_TYPE $threshold


## Add entries
{
	echo "GENE_ENV_TRAINING=$PHENO_DATA/geneEnvironmentTraining.csv"
	echo "GENE_ENV_TEST=$PHENO_DATA/geneEnvironmentTest.csv"
	echo "GENE_ENV_HOLDOUT=$PHENO_DATA/geneEnvironmentHoldout.csv"
	echo "MAIN_ENV_TRAINING=$PHENO_DATA/mainEnvironmentTraining.csv"
	echo "MAIN_ENV_TEST=$PHENO_DATA/mainEnvironmentTest.csv"
	echo "MAIN_ENV_HOLDOUT=$PHENO_DATA/mainEnvironmentHoldout.csv"
	
} >> "$CONFIG_FILE"

source $CONFIG_FILE

#exports done at top of script: PHENO_DATA, PHENO_PATH
bash "$SCRIPTS_DIR/run_plink_LD.sh" $PHENO

# ensure reduced feature file is present
#check if file exists post GxGxE filtering
for file in "$PHENO_DATA/scores/importantFeaturesPostShap.csv"; do
	if [[ -f "$file" ]]; then
		echo "✓ File exists: $file"
		REDUCED_FEATURE_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.csv"
	else
		echo "✗ File missing: $file"
	fi
done


# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--results_path "$RESULTS_PATH" \
--pheno_data "$PHENO_DATA" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--training_env_main_file "$MAIN_ENV_TRAINING" \
--test_env_main_file "$MAIN_ENV_TEST" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" \
--reduced_feature_file "$REDUCED_FEATURE_FILE"\
--epi_combo "$EPI_COMBO" \
--threshold $threshold


# Add entries
{
	echo "FEATURE_SCORES_FILE=$PHENO_DATA/scores/featureScoresReducedFinalModel.csv"
	echo "FINAL_MODEL_SCORES=$PHENO_DATA/scores/modelScoresReducedFinalModel.csv"
	echo "FINAL_MODEL_PROBABILITIES=$PHENO_DATA/scores/predictProbsReducedFinalModel.csv"
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
export SCORES_PATH="$PHENO_DATA/scores"
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH

#calculated PRS with beta coefs
# will run LD if not done previously and prune features in LD
# will filter GxGxE features post modelling to prune non-additive from separate genetic models
bash run_prs_calculations.sh $pheno


source $CONFIG_FILE

#ensure the McNemar test has been completed

#calculates prscr_mix trained using validation set and calculated for holdout set
export COVAR_FILE
python "${SCRIPTS_DIR}/prscr_modeling.py"

export RESULTS_PATH
python "${SCRIPTS_DIR}/clincal_measure_performance.py"

#exports done previously: PHENO_DATA, FEATURE_SCORES_FILE
#calculate unique features driving risk cohorts
export THRESHOLD=1.0
python "${SCRIPTS_DIR}/calculate_top_features_in_cohort.py"

conda deactivate

echo "TEST WORKFLOW COMPLETED!!"










    