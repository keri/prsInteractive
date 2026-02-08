#!/bin/bash


platform=${3:-"local"} 
#pheno=$1
PHENO='type2Diabetes'
threshold=${2:-1.99}
#EPI_COMBO=${4:-"sum"}
EPI_COMBO=${4:-"prod"}

# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi


if [ "${EPI_COMBO}" == "sum" ]; then
    COMBO_FOLDER='summedEpi'
else
    COMBO_FOLDER='productEpi'
fi

PHENO_DATA="${RESULTS_PATH}/$PHENO/${COMBO_FOLDER}"

source "$PHENO_DATA/pheno.config"

############## ACTIVATE ENVIRONMENT  ##############
# Check if environment exists, create if it doesn't
ENV_NAME="ukb_env"
#create conda shell to activate
if [ "$platform" != "hpc" ]; then
    # Initialize conda for bash (if not already done)
    if ! command -v conda &> /dev/null; then
        echo "ERROR: conda not found. Please install conda/miniconda first."
        exit 1
    fi
    
else
    module load Miniconda3/4.9.2
    source $(conda info --base)/etc/profile.d/conda.sh
    
fi

eval "$(conda shell.bash hook)"

if ! conda env list | grep -q "^${ENV_NAME}\s"; then
    echo "[CONDA] Creating conda environment from environment.yml..."
    conda env create -f "$PROJECT_ROOT/environment.yml"
    
else
    echo "[CONDA] Environment $ENV_NAME already exists"
    # Optionally update environment
    echo "[CONDA] Updating environment from environment.yml..."
    conda env update -f "$PRS_INTERACTIVE_HOME/environment.yml" --prune
fi

# Only activate conda environment if not running on HPC
if [ "$platform" != "hpc" ]; then
    echo "[SETUP] Local/standard platform detected - activating conda environment"
    
    # Activate the environment
    echo "[CONDA] Activating environment: $ENV_NAME"
    conda activate $ENV_NAME
    
    # Verify key tools are available
    echo "[CONDA] Verifying environment setup..."
    python --version
    R --version
    echo "[CONDA] Environment activated successfully"
else
    echo "[SETUP] HPC platform detected - conda environment created but not activated"
    echo "[SETUP] You can activate it later with: conda activate $ENV_NAME"
    echo "[SETUP] Or use system/module-provided tools as needed"
fi

required_files=(
    "$TRAINING_PATH"
    "$TEST_PATH"
    "$GENE_ENV_TRAINING"
    "$GENE_ENV_TEST"
    "$MAIN_ENV_TRAINING"
    "$MAIN_ENV_TEST"
    "$CLINICAL_ENV_TRAINING"
    "$CLINICAL_ENV_TEST"
    "$PHENO_DATA/scores/importantFeaturesPostShap.csv"
    "$PHENO_DATA/scores/cardioMetabolicimportantFeaturesPostShapFilteredZscore.csv"
    "${SCRIPTS_DIR}/glmPenalizedFinalModelling.R"
    
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file does not exist: $file"
        exit 1
    else
        echo "[DEBUG] Found: $file"
    fi
done
    
export PHENO_PATH=$PHENO_PATH
#export FEATURES_TO_FILTER_LD="$PHENO_DATA/scores/importantFeaturesPostShap.csv"
#check to see if LD has been done previously before association
if [ ! -f "$PHENO_DATA/scores/importantFeaturesPostShap.filteredLD.csv" ]; then
    echo "[LD FILTERING] LD is being done before modelling .."
    
    # Check for the actual PLINK output files
    if [ -f "$PHENO_PATH/finalModel.prune.in" ] || [ -f "$PHENO_PATH/finalModel.tags.list" ]; then
        echo "[LD FILTERING] Using existing LD files..."
        python "$SCRIPTS_DIR/filter_features_LD.py"
    else
        echo "[LD FILTERING] Creating new LD files..."
        python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"
        wait 
        
        plink --bfile "$PHENO_PATH/merged_allChromosomes" \
        --extract "$PHENO_PATH/finalModelLDSnps.txt" \
        --indep-pairwise 100kb 1 .6 \
        --r2 --show-tags all \
        --out "$PHENO_PATH/finalModel"
        
        python "$SCRIPTS_DIR/filter_features_LD.py"
    fi
    
    #check to see that reduced file script ran 
    if [ -f "$PHENO_DATA/scores/importantFeaturesPostShap.filteredLD.csv" ];then
        # Add entries
        {
            echo "REDUCED_FEATURE_FILE=$PHENO_DATA/scores/importantFeaturesPostShap.filteredLD.csv"
            
        } >> "${PHENO_DATA}/pheno.config"
        REDUCED_FEATURE_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.filteredLD.csv"
        
    else
        echo "[DEBUG] LD filtering didn't produce correct file"
        REDUCED_FEATURE_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.csv"
    fi
fi




REDUCED_FEATURE_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.filteredLD.csv"
### CHECK TO SEE IF LD HAS BEEN DONE

# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
--results_path "$RESULTS_PATH" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--pheno_data "$PHENO_DATA" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" \
--reduced_feature_file "$REDUCED_FEATURE_FILE" \
--epi_combo "$EPI_COMBO" \
--threshold $threshold \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--training_env_main_file "$MAIN_ENV_TRAINING" \
--test_env_main_file "$MAIN_ENV_TEST" \
--training_env_clinical_file "$CLINICAL_ENV_TRAINING" \
--test_env_clinical_file "$CLINICAL_ENV_TEST" 



# Add entries
{
    echo "FEATURE_SCORES_FILE=$PHENO_DATA/scores/featureScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_SCORES=$PHENO_DATA/scores/modelScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_PROBABILITIES=$PHENO_DATA/scores/predictProbsReducedFinalModel.csv"
} >> "${PHENO_DATA}/pheno.config"
    
filter non-additive gene-env features
python "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py" \
--config_file "$PHENO_DATA/pheno.config" \
--feature_scores_file "$PHENO_DATA/scores/featureScoresReducedFinalModel.csv"

#check to see if filtered file is present
if [ ! -f "$PHENO_DATA/scores/featureScoresReducedFinalModel.filtered.csv" ]; then
    echo "[G-E non-additive FILTERING] is not done"
    echo "run filter_non_additive_gene_env_features.py"
    exit 1
fi


python "$SCRIPTS_DIR/calculate_prs_post_modelling.py" \
--pheno_data $PHENO_DATA \
--test_file $TEST_PATH \
--holdout_file $HOLDOUT_PATH \
--covar_file $COVAR_FILE \
--hla_file $HLA_FILE \
--test_env_gen_file $GENE_ENV_TEST \
--holdout_env_gen_file $GENE_ENV_HOLDOUT \
--pheno $PHENO \
--feature_scores_file_filtered "$PHENO_DATA/scores/featureScoresReducedFinalModel.filtered.csv" \
--withdrawal_path $WITHDRAWAL_PATH \
--epi_combo $EPI_COMBO
bash run_prs_calculations.sh $PHENO $EPI_COMBO

python "$SCRIPTS_DIR/combine_prs" \
--pheno_data $PHENO_DATA

#calculate peformance of trained models and PRS calculations for main v other
#need PHENO_DATA exported or added as an --pheno_data argument
python "${SCRIPTS_DIR}/calculate_prs_stats.py" \
--pheno_data $PHENO_DATA

RScript "$SCRIPTS_DIR/prsAUCDelongStats.R" \
--pheno_data $PHENO_DATA


