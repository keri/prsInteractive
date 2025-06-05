#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

echo "=== Testing the EXACT mount commands from your Cromwell config ==="

# Test 1: Standard mounting strategy from your config
echo ""
echo "1. Testing standard mounting (as configured in your Cromwell config):"
docker run --rm \
    -v "$PROJECT_ROOT/data:/prsInteractive/data:ro" \
    -v "$PROJECT_ROOT/results:/prsInteractive/results" \
    -v "$PROJECT_ROOT/logs:/prsInteractive/logs" \
    -v "$PROJECT_ROOT/config:/prsInteractive/config:ro" \
    ubuntu:latest \
    /bin/bash -c "
        echo 'Standard mount test:';
        echo 'Contents of /prsInteractive/data:';
        ls -la /prsInteractive/data/ 2>/dev/null || echo 'FAILED: /prsInteractive/data not accessible';
        echo '';
        echo 'Checking for participant.csv:';
        if [ -f '/prsInteractive/data/participant.csv' ]; then
            echo 'SUCCESS: participant.csv found';
            ls -la /prsInteractive/data/participant.csv;
        else
            echo 'FAILED: participant.csv not found';
        fi
    "

echo ""
echo "2. Testing alternative mounting (whole project as host_prsInteractive):"
docker run --rm \
    -v "$PROJECT_ROOT:/host_prsInteractive:ro" \
    ubuntu:latest \
    /bin/bash -c "
        echo 'Alternative mount test:';
        echo 'Contents of /host_prsInteractive/data:';
        ls -la /host_prsInteractive/data/ 2>/dev/null || echo 'FAILED: /host_prsInteractive/data not accessible';
        echo '';
        echo 'Checking for participant.csv:';
        if [ -f '/host_prsInteractive/data/participant.csv' ]; then
            echo 'SUCCESS: participant.csv found via alternative mount';
            ls -la /host_prsInteractive/data/participant.csv;
        else
            echo 'FAILED: participant.csv not found via alternative mount';
        fi
    "

echo ""
echo "3. Testing with your actual Docker image (ukb-base:V1):"
if docker image inspect ukb-base:V1 >/dev/null 2>&1; then
    echo "Found ukb-base:V1 image, testing..."
    docker run --rm \
        -v "$PROJECT_ROOT/data:/prsInteractive/data:ro" \
        ukb-base:V1 \
        /bin/bash -c "
            echo 'Testing with ukb-base:V1:';
            echo 'Python version:';
            python --version 2>/dev/null || echo 'Python not found';
            echo '';
            echo 'Conda environment:';
            echo \$CONDA_DEFAULT_ENV;
            echo '';
            echo 'Can access data files:';
            ls -la /prsInteractive/data/participant.csv 2>/dev/null && echo 'SUCCESS: Files accessible in ukb-base:V1' || echo 'FAILED: Files not accessible in ukb-base:V1';
        "
else
    echo "ukb-base:V1 image not found - this might be part of the issue"
    echo "Available Docker images:"
    docker images | grep -E "(ukb|base)" || echo "No ukb-base images found"
fi

echo ""
echo "4. Testing the exact paths that will be used in the workflow:"
echo "Data path: $PROJECT_ROOT/data"
echo "Results path: $PROJECT_ROOT/results"
echo ""
echo "Files that should be accessible:"
for file in participant.csv participant_environment.csv covar.txt; do
    if [ -f "$PROJECT_ROOT/data/$file" ]; then
        echo "✅ $file exists on host"
    else
        echo "❌ $file missing on host"
    fi
done

echo ""
echo "5. Testing Cromwell execution directory access:"
mkdir -p "$PROJECT_ROOT/cromwell-executions/test"
echo "test file" > "$PROJECT_ROOT/cromwell-executions/test/testfile.txt"

docker run --rm \
    -v "$PROJECT_ROOT/cromwell-executions/test:/cromwell-executions" \
    -v "$PROJECT_ROOT/data:/prsInteractive/data:ro" \
    ubuntu:latest \
    /bin/bash -c "
        echo 'Testing cromwell execution directory:';
        ls -la /cromwell-executions/ || echo 'FAILED: cromwell-executions not accessible';
        echo '';
        echo 'Working directory simulation:';
        cd /cromwell-executions;
        pwd;
        echo 'Can access data from working directory:';
        ls -la /prsInteractive/data/participant.csv || echo 'FAILED: Cannot access data from working dir';
    "

# Cleanup
rm -f "$PROJECT_ROOT/cromwell-executions/test/testfile.txt"
rmdir "$PROJECT_ROOT/cromwell-executions/test" 2>/dev/null