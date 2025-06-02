#!/bin/bash
# Keri Multerer May 2025


# Set fixed project root directory
#PROJECT_ROOT="/nfs/scratch/projects/ukbiobank/prsInteractive"
#PROJECT_ROOT="$(dirname "$(pwd)")"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Export environment variable for Docker mounting
export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"


# Set base directories
#WORKFLOW_DIR="$PROJECT_ROOT/workflow"
#SCRIPTS_DIR="$PROJECT_ROOT/scripts"
#DATA_DIR="$PROJECT_ROOT/data"
#RESULTS_DIR="$PROJECT_ROOT/results"
#HPC_DIR="$PROJECT_ROOT/hpc"
#ENV_PATH="$PROJECT_ROOT/ukb_env"

# Create necessary directories
mkdir -p "$PROJECT_ROOT/results"
mkdir -p "$PROJECT_ROOT/logs"
mkdir -p "$PROJECT_ROOT/logs/call-logs"
mkdir -p "$PROJECT_ROOT/cromwell-executions"


# Create a .env file for future reference
cat > "$PROJECT_ROOT/.env" << EOF
# prsInteractive Pipeline Environment
PRS_INTERACTIVE_HOME=$PROJECT_ROOT
CROMWELL_ROOT=$PROJECT_ROOT/cromwell-executions
RESULTS_PATH=$PROJECT_ROOT/results
LOGS_DIR=$PROJECT_ROOT/logs
WORKFLOW_DIR="$PROJECT_ROOT/workflow"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_PATH="$PROJECT_ROOT/data"
HPC_DIR="$PROJECT_ROOT/hpc"
ENV_PATH="$PROJECT_ROOT/ukb_env"
EOF

# Set environment variable
#export DATA_PATH="$DATA_DIR"
#export RESULTS_PATH="$RESULTS_DIR"
#export ENV_PATH
#export SCRIPTS_DIR
#export HPC_DIR
#export PROJECT_ROOT


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "CONDA ENV BEING ACTIVATED ...: $ENV_PATH"
echo "[WORKFLOW] HPC_DIR IS SET TO: $HPC_DIR"


echo "Environment setup complete!"
echo ""
echo "To run the pipeline:"
echo "1. Source the environment: source .env"
echo "2. Run cromwell: java -jar \$CROMWELL_JAR run workflows/prsInteractivePipeline.DataProcessing.wdl -i workflows/pipelineInputs.json -o workflows/options.json --conf workflows/cromwell.conf"
echo ""
echo "Or use the run script: ./run_data_cleaning_workflow.sh"


