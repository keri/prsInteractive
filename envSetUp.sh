#!/bin/bash

#set up environment variables for use in workflow
platform=${5:-"local"}  # default to local
pheno=$1
icd10=$2
phenoStr=$3
n=$4 

set -e

# Generate a fixed configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# Create necessary pheno specific directories
echo "[SETUP] Creating directory structure..."
mkdir -p "$PROJECT_ROOT/results/$pheno"
mkdir -p "$PROJECT_ROOT/results/$pheno/scores"
mkdir -p "$PROJECT_ROOT/results/$pheno/models"
mkdir -p "$PROJECT_ROOT/results/$pheno/epiFiles"
mkdir -p "$PROJECT_ROOT/results/$pheno/epiFiles/preSummaryFiles"
mkdir -p "$PROJECT_ROOT/results/$pheno/figures"

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

# Create a .env file for future reference
cat > "$PROJECT_ROOT/results/$pheno/pheno.config" << EOF
# prsInteractive Pipeline Environment for $pheno
PHENO_PATH=$PROJECT_ROOT/results/$pheno
PHENO=$pheno
PHENO_STR="${phenoStr}"
ICD10=$icd10
N_CORES=$n
EPI_PATH=$PROJECT_ROOT/results/$pheno/epiFiles
PLATFORM=$platform
EOF

chmod 644 "$PROJECT_ROOT/results/$pheno/pheno.config"

echo "[WORKFLOW] DATA_PATH is set to: $PROJECT_ROOT/data"
echo "[WORKFLOW] RESULTS_PATH is set to: $PROJECT_ROOT/results"
echo "[WORKFLOW] Scripts directory: $PROJECT_ROOT/scripts"
echo "[WORKFLOW] PLATFORM: $platform"
echo "[WORKFLOW] CONDA ENV created: $ENV_NAME"
if [ "$platform" != "hpc" ]; then
    echo "[WORKFLOW] CONDA ENV activated: $ENV_NAME"
else
    echo "[WORKFLOW] Running on HPC - conda env created but not activated"
fi
echo "[WORKFLOW] HPC_DIR IS SET TO: $PROJECT_ROOT/hpc"
echo "[WORKFLOW] PHENOTYPE PATH is set to: $PROJECT_ROOT/results/$pheno"
echo "[PHENO] PHENO is set to: $pheno"
echo "[WORKFLOW] OUTPUT SCORES directory: $PROJECT_ROOT/results/$pheno/scores"
echo "[WORKFLOW] OUTPUT MODELS directory: $PROJECT_ROOT/results/$pheno/models"
echo "[WORKFLOW] OUTPUT EPI RESULTS directory: $PROJECT_ROOT/results/$pheno/epiFiles"
echo "[WORKFLOW] OUTPUT PHENO FIGURES directory: $PROJECT_ROOT/results/$pheno/figures"

DATA_PATH="$PROJECT_ROOT/data"
# Make sure data directory exists
if [ ! -d "$DATA_PATH" ]; then
    echo "WARNING: Data directory $DATA_PATH does not exist!"
    echo "Please create it and place your data files there: participant.csv participant_environment.csv covar.csv ukb_hla_v2.txt hla_participant.csv withdrawals.csv"
fi

# Check for required data files
echo "Checking for required data files..."
for file in "participant.csv" "participant_environment.csv" "ukb_hla_v2.txt" "hla_participant.csv" "withdrawals.csv"; do
    if [ ! -f "$DATA_PATH/$file" ]; then
        echo "WARNING: Required file $PROJECT_ROOT/data/$file not found!"
    else
        echo "✓ Found: $file"
    fi
done

echo "[SETUP] Environment setup complete!"
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
    