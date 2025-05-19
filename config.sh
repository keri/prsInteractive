#!/bin/bash
# Keri Multerer May 2025


# Set fixed project root directory
#PROJECT_ROOT="/nfs/scratch/projects/ukbiobank/prsInteractive"
PROJECT_ROOT="$(dirname "$(pwd)")"


# Set base directories
WORKFLOW_DIR="$PROJECT_ROOT/workflow"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
HPC_DIR="$PROJECT_ROOT/hpc"
ENV_PATH="$PROJECT_ROOT/ukb_env"

# Set environment variable
export DATA_PATH="$DATA_DIR"
export RESULTS_PATH="$RESULTS_DIR"
export ENV_PATH
export SCRIPTS_DIR
export HPC_DIR
export PROJECT_ROOT


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "CONDA ENV BEING ACTIVATED ...: $ENV_PATH"
echo "[WORKFLOW] HPC_DIR IS SET TO: $HPC_DIR"



