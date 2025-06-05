#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Testing Cromwell Configuration ==="

# Source the environment
if [ -f "$PROJECT_ROOT/.env" ]; then
    source "$PROJECT_ROOT/.env"
    echo "Environment loaded:"
    echo "CROMWELL_CONFIG: $CROMWELL_CONFIG"
    echo "PRS_INTERACTIVE_HOME: $PRS_INTERACTIVE_HOME"
else
    echo "❌ No .env file found"
    exit 1
fi

if [ -z "$CROMWELL_CONFIG" ]; then
    echo "❌ CROMWELL_CONFIG is empty"
    exit 1
fi

if [ ! -f "$CROMWELL_CONFIG" ]; then
    echo "❌ Config file doesn't exist: $CROMWELL_CONFIG"
    exit 1
fi

echo ""
echo "=== Testing Docker command from config ==="

# Extract the docker command from the config
echo "Config file contents:"
cat "$CROMWELL_CONFIG"

echo ""
echo "=== Manual test of the Docker mount command ==="

# Test the exact command that Cromwell should be using
ABS_PROJECT_ROOT=$(realpath "$PROJECT_ROOT")
echo "Absolute project root: $ABS_PROJECT_ROOT"

# Test the mount command manually
docker run --rm \
-v "$ABS_PROJECT_ROOT/data:/prsInteractive/data:ro" \
-v "$ABS_PROJECT_ROOT/results:/prsInteractive/results" \
-v "$ABS_PROJECT_ROOT/logs:/prsInteractive/logs" \
-v "$ABS_PROJECT_ROOT/config:/prsInteractive/config:ro" \
ukb-base:V1 \
/bin/bash -c "
        echo 'Testing mount in ukb-base:V1 container:';
        echo 'Contents of /prsInteractive/data:';
        ls -la /prsInteractive/data/ || echo '❌ /prsInteractive/data not accessible';
        echo '';
        echo 'Testing participant.csv access:';
        python3 -c \"
import pandas as pd
try:
    df = pd.read_csv('/prsInteractive/data/participant.csv', nrows=5)
    print('✅ SUCCESS: Python can read participant.csv')
    print(f'Columns: {list(df.columns)[:5]}')
except Exception as e:
    print(f'❌ FAILED: {e}')
\"
    "

echo ""
echo "=== If that worked, the issue is in how Cromwell is calling Docker ==="
echo "Let's check if we need to update the Cromwell command itself..."

# Show how to run Cromwell with the correct config
echo ""
echo "=== Correct Cromwell execution command ==="
echo "To run Cromwell with the proper configuration:"
echo ""
echo "java -jar cromwell.jar run \\"
echo "    workflows/prsInteractivePipeline.wdl \\"
echo "    --inputs workflows/pipelineInputs.json \\"
echo "    --options $CROMWELL_CONFIG"
echo ""
echo "OR if using a separate options file:"
echo ""
echo "java -jar cromwell.jar run \\"
echo "    workflows/prsInteractivePipeline.wdl \\"
echo "    --inputs workflows/pipelineInputs.json \\"
echo "    --options cromwell-options.json"

