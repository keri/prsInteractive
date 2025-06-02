#!/bin/bash

# Portable run script for prsInteractive pipeline

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Source environment if it exists
if [ -f ".env" ]; then
    source .env
    echo "Loaded environment from .env"
else
    echo "No .env file found, running setup..."
    echo "create environment by running envSetUp.sh <pheno> <icd10> <pheno string> <n cores for epistatic analysis> "
    exit 1
fi


# Set default cromwell jar location if not specified
CROMWELL_JAR="${CROMWELL_JAR:-$HOME/cromwell/cromwell-90.jar}"

# Check if cromwell jar exists
if [ ! -f "$CROMWELL_JAR" ]; then
    echo "Error: Cromwell JAR not found at $CROMWELL_JAR"
    echo "Please set CROMWELL_JAR environment variable or place cromwell jar at the default location"
    exit 1
fi

# Export the project root for Docker mounting
export PRS_INTERACTIVE_HOME="$SCRIPT_DIR"

echo "Running prsInteractive pipeline..."
echo "Project root: $PRS_INTERACTIVE_HOME"
echo "Cromwell JAR: $CROMWELL_JAR"

# Run the workflow
java -Dconfig.file=workflows/cromwell.conf -jar "$CROMWELL_JAR" run \
    workflows/prsInteractivePipeline.DataProcessing.wdl \
    -i workflows/pipelineInputs.json \
    -o workflows/options.json

echo "Pipeline execution completed!"