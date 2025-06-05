#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Debugging Cromwell Configuration Loading ==="

# Load environment
source .env

echo "Environment variables:"
echo "  CROMWELL_CONFIG: $CROMWELL_CONFIG"
echo "  PLATFORM: $PLATFORM"

echo ""
echo "=== Checking configuration file content ==="
if [ -f "$CROMWELL_CONFIG" ]; then
    echo "✅ Config file exists: $CROMWELL_CONFIG"
    echo "Content:"
    cat "$CROMWELL_CONFIG"
else
    echo "❌ Config file missing: $CROMWELL_CONFIG"
fi

echo ""
echo "=== Checking default.config (referenced by local.config) ==="
DEFAULT_CONFIG="$(dirname "$CROMWELL_CONFIG")/default.config"
if [ -f "$DEFAULT_CONFIG" ]; then
    echo "✅ Default config exists: $DEFAULT_CONFIG"
    echo "Content (first 50 lines):"
    head -50 "$DEFAULT_CONFIG"
else
    echo "❌ Default config missing: $DEFAULT_CONFIG"
fi

echo ""
echo "=== Testing the Docker command from config manually ==="
ABS_PROJECT_ROOT=$(realpath "$PROJECT_ROOT")

# Extract the docker command that should be used
echo "Expected Docker command based on config:"
echo "docker run --rm -v \${cwd}:\${docker_cwd} -v $ABS_PROJECT_ROOT/data:/prsInteractive/data:ro -v $ABS_PROJECT_ROOT/results:/prsInteractive/results ukb-base:V1 /bin/bash \${script}"

echo ""
echo "=== Testing actual volume mount manually ==="
docker run --rm \
-v "$ABS_PROJECT_ROOT/data:/prsInteractive/data:ro" \
-v "$ABS_PROJECT_ROOT/results:/prsInteractive/results" \
ukb-base:V1 \
/bin/bash -c "echo 'Manual test:'; ls -la /prsInteractive/data/participant.csv || echo 'FAILED: File not accessible'"

echo ""
echo "=== Checking the actual Docker command Cromwell generated ==="
LATEST_EXEC=$(find "$PROJECT_ROOT/cromwell-executions" -name "execution" -type d | tail -1)
if [ -n "$LATEST_EXEC" ] && [ -d "$LATEST_EXEC" ]; then
    echo "Latest execution directory: $LATEST_EXEC"
    
    if [ -f "$LATEST_EXEC/script" ]; then
        echo ""
        echo "Generated script content:"
        cat "$LATEST_EXEC/script"
    fi
    
    # Look for Docker command in Cromwell logs
    echo ""
    echo "Looking for Docker command in execution artifacts..."
    find "$LATEST_EXEC/.." -name "*.log" -exec grep -l "docker run" {} \; 2>/dev/null | head -3 | while read logfile; do
        echo "Docker command from $logfile:"
        grep "docker run" "$logfile" || echo "No docker run command found"
    done
fi

echo ""
echo "=== Cromwell version and configuration test ==="
# Test if Cromwell recognizes our config
CROMWELL_JAR="${CROMWELL_JAR:-$HOME/cromwell/cromwell-90.jar}"
if [ -f "$CROMWELL_JAR" ]; then
    echo "Testing Cromwell config parsing..."
    echo "java -Dconfig.file=$CROMWELL_CONFIG -jar $CROMWELL_JAR --help" 
    timeout 10s java -Dconfig.file="$CROMWELL_CONFIG" -jar "$CROMWELL_JAR" --help 2>&1 | head -5 || echo "Config test timed out or failed"
fi

echo ""
echo "=== Diagnosis ==="
echo "If manual Docker test works but Cromwell doesn't:"
echo "1. Cromwell may not be loading the config file properly"
echo "2. The -Dconfig.file syntax might be wrong for your Cromwell version"
echo "3. There might be a syntax error in the config file"
echo "4. Cromwell might be using a default backend instead of 'Local'"

