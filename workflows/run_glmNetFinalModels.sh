#!/bin/bash


platform=${3:-"local"} 
pheno=$1
#pheno='type2Diabetes'
threshold=${2:-1.99}
EPI_COMBO=${4:-"sum"}

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
    conda env update -f "$PROJECT_ROOT/environment.yml" --prune
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



# ensure reduced feature file is present
for file in "$PHENO_DATA/scores/importantFeaturesPostShap.csv"; do
    if [[ -f "$file" ]]; then
        echo "✓ File exists: $file"
    else
        echo "✗ File missing: $file"
    fi
done

#check to see if LD has been done previously before association
if [ ! -f "$PHENO_DATA/finalModel.ld" ];then
    export PHENO_DATA
    export PHENO_PATH
    bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
fi

REDUCED_FEATURE_FILE="$PHENO_DATA/scores/importantFeaturesPostShap.csv"

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

# Pass to R script
Rscript "$SCRIPTS_DIR/glmPenalizedFinalModelling.R" \
--results_path "$RESULTS_PATH" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--pheno_data "$PHENO_DATA" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--training_env_main_file "$MAIN_ENV_TRAINING" \
--test_env_main_file "$MAIN_ENV_TEST" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" \
--reduced_feature_file "$REDUCED_FEATURE_FILE" \
--epi_combo "$EPI_COMBO"\
--threshold $threshold


# Add entries
{
    echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_SCORES=$PHENO_PATH/scores/modelScoresReducedFinalModel.csv"
    echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scores/predictProbsReducedFinalModel.csv"
} >> "${PHENO_DATA}/pheno.config"
    
bash run_prs_calculations.sh $pheno

#calculate peformance of trained models and PRS calculations for main v other
#need PHENO_PATH exported or added as an --pheno_path argument
python "${SCRIPTS_DIR}/calculate_prs_stats.py"