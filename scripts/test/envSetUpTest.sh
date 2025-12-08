#!/bin/bash

#set up environment variables for use in workflow
platform=${6:-"local"}  # default to local
pheno=$1
icd10=$2
phenoStr=$3
n=$4 
EPI_COMBO=${5:-"sum"}


set -e

# Generate a fixed configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Export environment variable for Docker mounting
export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"

echo "[SETUP] Setting up environment for phenotype: $pheno"
echo "[SETUP] Platform: $platform"

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



# Create a .env file for future reference
cat > "$PROJECT_ROOT/env.config" << EOF
# prsInteractive Pipeline Environment
PRS_INTERACTIVE_HOME=$PROJECT_ROOT
RESULTS_PATH=$PROJECT_ROOT/results
WORKFLOW_DIR=$PROJECT_ROOT/workflow
SCRIPTS_DIR=$PROJECT_ROOT/scripts
DATA_PATH=$PROJECT_ROOT/data
HPC_DIR=$PROJECT_ROOT/hpc
ENV_PATH=$PROJECT_ROOT/ukb_env
CONDA_ENV_NAME=$ENV_NAME
PLATFORM=$platform
WITHDRAWAL_PATH=$PROJECT_ROOT/data/withdrawals.csv
EOF

chmod 644 "$PROJECT_ROOT/env.config"

source "$PROJECT_ROOT/env.config"


#### CREATE PHENO VARIABLES  ################

# Create folders to store results based on how epi is combined
# Set PHENO_PATH based on EPI_COMBO value
if [ "${EPI_COMBO}" == "sum" ]; then
    COMBO_FOLDER='summedEpi'
else
    COMBO_FOLDER='productEpi'
fi
PHENO_DATA="${RESULTS_PATH}/$pheno/${COMBO_FOLDER}"
echo "[SETUP] Creating directory structure..."

mkdir -p "$PROJECT_ROOT/results/$pheno"
mkdir -p "$PHENO_DATA"
mkdir -p "$PHENO_DATA/scores"
mkdir -p "$PHENO_DATA/models"
mkdir -p "$PHENO_DATA/figures"

#make folders present for any epi combo
mkdir -p "$PROJECT_ROOT/results/$pheno/epiFiles"
mkdir -p "$PROJECT_ROOT/results/$pheno/epiFiles/preSummaryFiles"

# Create necessary pheno specific directories
echo "[SETUP] Creating directory structure..."

mkdir -p "$PROJECT_ROOT/testData"
mkdir -p "$PROJECT_ROOT/testData/variant_calls"
mkdir -p "$PROJECT_ROOT/results/ParticipantTestData"

#provide permissions for folder
chmod +x "$PROJECT_ROOT/testData/"
chmod +x "$PROJECT_ROOT/testData/variant_calls/"
chmod +x "$PHENO_DATA/scores/"
chmod +x "$PHENO_DATA/models/"
chmod +x "$PHENO_DATA/figures/"
chmod +x "$PROJECT_ROOT/results/$pheno/epiFiles/"
chmod +x "$PROJECT_ROOT/results/$pheno/epiFiles/preSummaryFiles/"
chmod +x "$PROJECT_ROOT/results/ParticipantTestData/"

# Create a .env file for future reference
cat > "$PHENO_DATA/pheno.config" << EOF
# prsInteractive Pipeline Environment for $pheno
PHENO_PATH="$PROJECT_ROOT/results/$pheno"
RESULTS_PATH="$PROJECT_ROOT/results/ParticipantTestData"
PHENO_DATA=$PHENO_DATA
PHENO=$pheno
PHENO_STR="${phenoStr}"
ICD10=$icd10
N_CORES=$n
EPI_PATH=$PROJECT_ROOT/results/$pheno/epiFiles
PLATFORM=$platform
WITHDRAWAL_PATH="$PROJECT_ROOT/testData/withdrawals.csv"
DATA_PATH="$PROJECT_ROOT/testData"
EOF

chmod 644 "$PHENO_DATA/pheno.config"

source "$PHENO_DATA/pheno.config"

echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "[WORKFLOW] PLATFORM: $platform"
echo "[WORKFLOW] CONDA ENV created: $ENV_NAME"
if [ "$platform" != "hpc" ]; then
    echo "[WORKFLOW] CONDA ENV activated: $ENV_NAME"
else
    echo "[WORKFLOW] Running on HPC - conda env created but not activated"
fi
echo "[WORKFLOW] HPC_DIR IS SET TO: $HPC_DIR"
echo "[WORKFLOW] PHENOTYPE PATH is set to: $PHENO_PATH"
echo "[PHENO] PHENO is set to: $pheno"
echo "[WORKFLOW] OUTPUT SCORES directory: $PHENO_DATA/scores"
echo "[WORKFLOW] OUTPUT MODELS directory: $PHENO_DATA/models"
echo "[WORKFLOW] OUTPUT EPI RESULTS directory: $PHENO_PATH/epiFiles"
echo "[WORKFLOW] OUTPUT PHENO FIGURES directory: $PHENO_DATA/figures"


# Make sure data directory exists
if [ ! -d "$DATA_PATH" ]; then
    echo "WARNING: Data directory $DATA_PATH does not exist!"
fi


#############################  CREATE TEST DATA ######################################

export DATA_PATH
export SCRIPTS_DIR
#create genotype data in .bed format
bash "$SCRIPTS_DIR/test/create_simulation_genotype_data.sh"


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/test/create_simulated_participant_covariate_data.py"



export PHENO=$pheno 
export PHENO_PATH
export PHENO_STR
export icd10
# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"

export RESULTS_PATH
# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Add entries
{
    echo "HLA_FILE=$RESULTS_PATH/participant_hla.csv"
    echo "COVAR_FILE=$RESULTS_PATH/covar.csv" 
    echo "ENV_FILE=$RESULTS_PATH/participant_environment.csv"
} >> "$PHENO_DATA/pheno.config"

# Simple verification
echo "Config file updated. Contents:"
cat "$PHENO_DATA/pheno.config"


# Check for required cleaned data files
echo "Checking for required data files..."
for file in  "participant_environment.csv" "covar.csv" "participant_hla.csv"; do
    if [ ! -f "$RESULTS_PATH/$file" ]; then
        echo "WARNING: Required file $RESULTS_PATH/$file not found!"
    else
        echo "✓ Found: $file"
    fi
done


# Run the variant call cleaning
# test script to handle special SNPs in chr 6
bash "$SCRIPTS_DIR/test/plink_clean_variant_calls_test.sh"

echo "[SETUP] Environment setup complete!"
echo "[SETUP] Test data created!"

echo ""
if [ "$platform" != "hpc" ]; then
    echo "To run workflow steps, use:"
    echo "  source env.config && conda activate $ENV_NAME"
    echo "  # then run your scripts"
else
    echo "HPC platform detected. To run workflow steps, choose one option:"
    echo "  Option 1 - Use conda environment:"
    echo "    source env.config && conda activate $ENV_NAME"
    echo "    # then run your scripts"
    echo "  Option 2 - Use system/module tools:"
    echo "    source env.config"
    echo "    # Load required modules (e.g., module load python R)"
    echo "    # then run your scripts"
fi

#bash "$PROJECT_ROOT/update_pipeline_inputs.sh" "$pheno" "$icd10" "$phenoStr" $n "$platform"
    