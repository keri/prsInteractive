#!/bin/bash

# Portable run script for prsInteractive pipeline
# Compatible with Cromwell 90

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== prsInteractive Pipeline Runner ==="
echo "Script directory: $SCRIPT_DIR"

# Source environment if it exists
if [ -f ".env" ]; then
    source .env
    echo "✅ Loaded environment from .env"
    
    # Export the variables you need in Docker
    export PRS_INTERACTIVE_HOME
    export RESULTS_PATH
    export SCRIPTS_DIR
    export DATA_PATH

    
    echo "Environment variables:"
    echo "  PRS_INTERACTIVE_HOME: $PRS_INTERACTIVE_HOME"
    echo "  RESULTS PATH : $RESULTS_PATH"

else
    echo "❌ No .env file found, running setup..."
    echo "Create environment by running:"
    echo "  ./envSetUp.sh <pheno> <icd10> <pheno_string> <n_cores> <platform>"
    echo "Example:"
    echo "  ./envSetUp.sh type2Diabetes E11 'type 2 diabetes' 40 local"
    exit 1
fi


required

#for file in "${required_files[@]}"; do
#   if [ -f "$file" ]; then
#       echo "✅ Found: $file"
#   else
#       echo "❌ Missing: $file"
#       exit 1
#   fi
#done

echo ""
echo "=== Pipeline Configuration ==="
echo "Project root: $PRS_INTERACTIVE_HOME"
echo "Cromwell JAR: $CROMWELL_JAR"
echo "Platform: $PLATFORM"
echo "Config file: $CROMWELL_CONFIG"
echo "WDL file: workflows/prsInteractivePipeline.wdl"
echo "Inputs: workflows/pipelineInputs.json"
echo "Options: workflows/options.json"


exit_code=$?

echo ""
echo "=== Pipeline Execution Completed ==="
echo "$(date): Pipeline finished with exit code: $exit_code"

if [ $exit_code -eq 0 ]; then
    echo "✅ SUCCESS: Pipeline completed successfully!"
    echo ""
    echo "Results should be available in:"
    echo "  $RESULTS_PATH"
    echo ""

else
    echo "❌ FAILED: Pipeline failed with exit code $exit_code"
    echo ""

fi

exit $exit_code