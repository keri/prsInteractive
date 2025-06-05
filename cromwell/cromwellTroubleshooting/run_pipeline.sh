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
    export PLATFORM
    export CROMWELL_CONFIG
    
    echo "Environment variables:"
    echo "  PRS_INTERACTIVE_HOME: $PRS_INTERACTIVE_HOME"
    echo "  PLATFORM: $PLATFORM"
    echo "  CROMWELL_CONFIG: $CROMWELL_CONFIG"
else
    echo "❌ No .env file found, running setup..."
    echo "Create environment by running:"
    echo "  ./envSetUp.sh <pheno> <icd10> <pheno_string> <n_cores> <platform>"
    echo "Example:"
    echo "  ./envSetUp.sh type2Diabetes E11 'type 2 diabetes' 40 local"
    exit 1
fi

# Set cromwell jar location - check multiple possibilities
CROMWELL_JAR_OPTIONS=(
    "${CROMWELL_JAR:-}"                           # User-specified
    "$HOME/cromwell/cromwell-90.jar"               # Your current version
    "$SCRIPT_DIR/cromwell.jar"                     # Local to project
    "$HOME/cromwell/cromwell-85.jar"               # Older version
    "$HOME/Downloads/cromwell-90.jar"              # Download location
)

CROMWELL_JAR=""
for jar in "${CROMWELL_JAR_OPTIONS[@]}"; do
    if [ -n "$jar" ] && [ -f "$jar" ]; then
        CROMWELL_JAR="$jar"
        break
    fi
done

# Check if we found cromwell
if [ -z "$CROMWELL_JAR" ] || [ ! -f "$CROMWELL_JAR" ]; then
    echo "❌ Cromwell JAR not found. Please ensure cromwell is available."
    echo "Checked locations:"
    for jar in "${CROMWELL_JAR_OPTIONS[@]}"; do
        echo "  - $jar"
    done
    exit 1
fi

# Detect Cromwell version for syntax compatibility
echo ""
echo "=== Detecting Cromwell Version ==="
CROMWELL_VERSION=$(java -jar "$CROMWELL_JAR" --version 2>/dev/null | head -1 || echo "unknown")
echo "Cromwell version: $CROMWELL_VERSION"

# Verify required files exist
echo ""
echo "=== Verifying required files ==="
required_files=(
    "workflows/prsInteractivePipeline.wdl"
    "workflows/pipelineInputs.json"
    "workflows/options.json"
    "$CROMWELL_CONFIG"
)

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ Found: $file"
    else
        echo "❌ Missing: $file"
        exit 1
    fi
done

echo ""
echo "=== Pipeline Configuration ==="
echo "Project root: $PRS_INTERACTIVE_HOME"
echo "Cromwell JAR: $CROMWELL_JAR"
echo "Platform: $PLATFORM"
echo "Config file: $CROMWELL_CONFIG"
echo "WDL file: workflows/prsInteractivePipeline.wdl"
echo "Inputs: workflows/pipelineInputs.json"
echo "Options: workflows/options.json"

# Build the Cromwell command based on version
# Cromwell 90+ uses -Dconfig.file syntax
CONFIG_OPTION="-Dconfig.file=$CROMWELL_CONFIG"

# Show the command that will be executed
echo ""
echo "=== Cromwell Command (v90+ syntax) ==="
echo "java -Xmx8g $CONFIG_OPTION -jar $CROMWELL_JAR run \\"
echo "    workflows/prsInteractivePipeline.wdl \\"
echo "    -i workflows/pipelineInputs.json \\"
echo "    -o workflows/options.json"

echo ""
echo "=== Starting Pipeline Execution ==="
echo "$(date): Starting prsInteractive pipeline..."

# Run the workflow with Cromwell 90 syntax
java -Xmx8g "$CONFIG_OPTION" -jar "$CROMWELL_JAR" run \
    workflows/prsInteractivePipeline.DataProcessing.wdl \
    -i workflows/pipelineInputs.json \
    -o workflows/options.json

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
    echo "Logs available in:"
    echo "  ./logs/"
else
    echo "❌ FAILED: Pipeline failed with exit code $exit_code"
    echo ""
    echo "Check the logs for details:"
    echo "  ./logs/"
    echo "  ./cromwell-executions/"
    
    # Show recent error logs if available
    echo ""
    echo "=== Recent Error Logs ==="
    if [ -d "./cromwell-executions" ]; then
        # Find the most recent stderr file
        LATEST_ERROR=$(find ./cromwell-executions -name "stderr" -type f -exec stat -f "%m %N" {} \; 2>/dev/null | sort -nr | head -1 | cut -d' ' -f2-)
        if [ -n "$LATEST_ERROR" ]; then
            echo "Latest error from: $LATEST_ERROR"
            echo "Last 10 lines:"
            tail -10 "$LATEST_ERROR"
        fi
    fi
fi

exit $exit_code