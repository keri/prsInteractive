#!/bin/bash

#set up environment variables for use in workflow

pheno=$1
icd10=$2
phenoStr=$3
n=$4 

#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"
#n=40

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "project root is : ${PROJECT_ROOT}"

# Export environment variable for Docker mounting
export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"

# Set base directories
#make a folder inside root data folder for each phenotype

# Create necessary directories
mkdir -p "$PROJECT_ROOT/results/$pheno"
mkdir -p "$PROJECT_ROOT/logs"
mkdir -p "$PROJECT_ROOT/logs/call-logs"
mkdir -p "$PROJECT_ROOT/cromwell-executions"

PHENO_DIR="$PROJECT_ROOT/results/$pheno"


echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $pheno"
echo "PHENOTYPE STRING USED IN ANALYSIS ...: $phenoStr"
echo "ICD10 USED IN ANALYSIS ...: $icd10"


#check to see if PHENO_PATH is already present
#create a phenotype env for later use
# Create a .env file for future reference
cat > "${PHENO_DIR}/pheno_config.sh" << EOF
#pheno environment variables
PHENO="${pheno}"
PHENO_PATH=$PHENO_DIR
PHENO_STR="${phenoStr}"
ICD10=$icd10
N=$n
EOF

# Create a .env file for future reference
cat > "$PROJECT_ROOT/.env" << EOF
# prsInteractive Pipeline Environment
PRS_INTERACTIVE_HOME=$PROJECT_ROOT
CROMWELL_ROOT=$PROJECT_ROOT/cromwell-executions
RESULTS_PATH=$PROJECT_ROOT/results
LOGS_DIR=$PROJECT_ROOT/logs
WORKFLOW_DIR=$PROJECT_ROOT/workflow
SCRIPTS_DIR=$PROJECT_ROOT/scripts
DATA_PATH=$PROJECT_ROOT/data
HPC_DIR=$PROJECT_ROOT/hpc
ENV_PATH=$PROJECT_ROOT/ukb_env
EOF


echo "[WORKFLOW] DATA_PATH is set to: $PROJECT_ROOT/data"
echo "[WORKFLOW] RESULTS_PATH is set to: $PROJECT_ROOT/results"
echo "[WORKFLOW] Scripts directory: $PROJECT_ROOT/scripts"
echo "CONDA ENV BEING ACTIVATED ...: $PROJECT_ROOT/ukb_en"
echo "[WORKFLOW] HPC_DIR IS SET TO: $PROJECT_ROOT/hpc"



