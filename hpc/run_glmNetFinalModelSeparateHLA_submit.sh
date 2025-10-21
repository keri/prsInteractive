#!/bin/bash

#
#SBATCH --job-name=glmFinalModel
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_glmFinalModel.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_glmFinalModel.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=70
#SBATCH --mem=800G
#SBATCH --time=90:00:00
#


#load required models
module load GCC/11.2.0  
module load OpenMPI/4.1.1
module load R/4.2.0

module load plink/1.90

pheno=$1
# Source config
source ../env.config # because you're in prsInteractive/hpc


#check that a results folder for phenotype exists
if [ ! -d "${RESULTS_PATH}/$pheno" ]; then
    echo "Folder '${RESULTS_PATH}/$pheno' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${RESULTS_PATH}/$pheno/pheno.config"
fi

export RESULTS_PATH=$RESULTS_PATH
export PHENO_PATH=$PHENO_PATH
export DATA_PATH=$DATA_PATH
export HLA_FILE="$RESULTS_PATH/participant_hla.csv"
export COVAR_FILE="$RESULTS_PATH/covar.csv"
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export TRAINING_ENV_GEN_FILE=$GENE_ENV_TRAINING
export TEST_ENV_GEN_FILE=$GENE_ENV_TEST

### CHECK TO SEE IF LD HAS BEEN DONE

if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
    echo "LD not completed for importantFeaturesPostShap.csv"
    #run LD script
    export PRE_POST_ASSOCIATION='pre'
    bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
fi


#check to see LD was completed
if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
    echo "LD analysis failed for importantFeaturesPostShap.csv"
    FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesPostShap.csv"
    echo "reverting to original importantFeaturesPostShap.csv"
else
    echo "LD analysis completed for importantFeaturesPostShap.csv"
    FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv"
fi

#bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
echo "Starting R script..."

# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModellingSeparateHLA.R" \
--results_path "$RESULTS_PATH" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--pheno_path "$PHENO_PATH" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--feature_model_file "$FEATURES_FOR_MODEL_FILE"

R_EXIT_CODE=$?

if [ $R_EXIT_CODE -ne 0 ]; then
    echo "ERROR: R script failed with exit code $R_EXIT_CODE"
    exit 1
fi

echo "R script completed successfully. Checking outputs..."

# Check if R script outputs exist before updating config
missing_files=()

if [ ! -f "$PHENO_PATH/scores/featureScoresReducedFinalModelSeparateHLA.csv" ]; then
    missing_files+=("featureScoresReducedFinalModelSeparateHLA.csv")
fi

if [ ! -f "$PHENO_PATH/scores/modelScoresReducedFinalModelSeparateHLA.csv" ]; then
    missing_files+=("modelScoresReducedFinalModelSeparateHLA.csv")
fi

if [ ! -f "$PHENO_PATH/scores/predictProbsReducedFinalModelSeparateHLA.csv" ]; then
    missing_files+=("predictProbsReducedFinalModelSeparateHLA.csv")
fi

# Only update config if all files exist
if [ ${#missing_files[@]} -eq 0 ]; then
    echo "All R script outputs found. Updating pheno.config..."
    {
        echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModelSeparateHLA.csv"
        echo "FINAL_MODEL_SCORES=$PHENO_PATH/scores/modelScoresReducedFinalModelSeparateHLA.csv"
        echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scores/predictProbsReducedFinalModelSeparateHLA.csv"
    } >> "${RESULTS_PATH}/$pheno/pheno.config"
else
    echo "ERROR: R script failed to generate the following files:"
    printf '%s\n' "${missing_files[@]}"
    echo "pheno.config not updated"
    exit 1
fi

#sbatch "run_prs_calculations_submit.sh" $pheno

