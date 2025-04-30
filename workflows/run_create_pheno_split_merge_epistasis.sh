#!/bin/bash

pheno=$1

# Set base directories
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$WORKFLOW_DIR")"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="$PROJECT_ROOT/data/"
PHENO_DIR="$PROJECT_ROOT/data/$pheno/"

# Set environment variable
export DATA_PATH="$DATA_DIR"

echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"


