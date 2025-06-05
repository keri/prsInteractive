#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# Load environment
source .env

CROMWELL_JAR="${CROMWELL_JAR:-$HOME/cromwell/cromwell-90.jar}"

echo "=== Running with Docker-Forced Configuration ==="
echo "Config: $PROJECT_ROOT/config/simple-docker.config"
echo "Command:"
echo "java -Dconfig.file=$PROJECT_ROOT/config/simple-docker.config -jar $CROMWELL_JAR run workflows/prsInteractivePipeline.DataProcessing.wdl -i workflows/pipelineInputs.json -o workflows/options.json"
echo ""

java -Dconfig.file="$PROJECT_ROOT/config/simple-docker.config" -jar "$CROMWELL_JAR" run \
    workflows/prsInteractivePipeline.DataProcessing.wdl \
    -i workflows/pipelineInputs.json \
    -o workflows/options.json
