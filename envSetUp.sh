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


# Create necessary pheno specific directories
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
WITHDRAWN_PATH=$PROJECT_ROOT/data/withdrawals.csv
EOF

chmod 777 "$PROJECT_ROOT/env.config"


# Create a .env file for future reference
cat > "$PROJECT_ROOT/results/$pheno/pheno.config" << EOF
# prsInteractive Pipeline Environment for $pheno
PHENO_PATH=$PROJECT_ROOT/results/$pheno
PHENO=$pheno
PHENO_STR="${phenoStr}"
ICD10=$icd10
N_CORES=$n
EPI_PATH=$PROJECT_ROOT/results/$pheno/epiFiles
EOF

chmod 777 "$PROJECT_ROOT/results/$pheno/pheno.config"

echo "[WORKFLOW] DATA_PATH is set to: $PROJECT_ROOT/data"
echo "[WORKFLOW] RESULTS_PATH is set to: $PROJECT_ROOT/results"
echo "[WORKFLOW] Scripts directory: $PROJECT_ROOT/scripts"
echo "CONDA ENV BEING ACTIVATED ...: $PROJECT_ROOT/ukb_env"
echo "[WORKFLOW] HPC_DIR IS SET TO: $PROJECT_ROOT/hpc"


echo "[WORKFLOW] PHENOTYPE PATH is set to: $PROJECT_ROOT/results/$pheno"
echo "[PHENO] PHENO is set to: $pheno"
echo "[WORKFLOW] OUTPUT SCORES directory: $PROJECT_ROOT/results/$pheno/scores"
echo "[WORKFLOW] OUTPUT MODELS directory:: $PROJECT_ROOT/results/$pheno/models"
echo "[WORKFLOW] OUTPUT EPI RESULTS directory:: $PROJECT_ROOT/results/$pheno/epiFiles"
echo "[WORKFLOW] OUTPUT PHENO FIGURES directory:: $PROJECT_ROOT/results/$pheno/figures"


DATA_PATH="$PROJECT_ROOT/data"
# Make sure data directory exists
if [ ! -d "$DATA_PATH" ]; then
    echo "WARNING: Data directory $DATA_PATH does not exist!"
    echo "Please create it and place your data files there: participant.csv participant_environment.csv covar.csv ukb_hla_v2.txt hla_participant.csv withdrawals.csv"
fi
# Check for required data files
echo "Checking for required data files..."
for file in "participant.csv" "participant_environment.csv" "covar.csv" "ukb_hla_v2.txt" "hla_participant.csv" "withdrawals.csv"; do
    if [ ! -f "$DATA_PATH/$file" ]; then
        echo "WARNING: Required file $PROJECT_ROOT/data/$file not found!"
    else
        echo "✓ Found: $file"
    fi
done
    



#bash "$PROJECT_ROOT/update_pipeline_inputs.sh" "$pheno" "$icd10" "$phenoStr" $n "$platform"


