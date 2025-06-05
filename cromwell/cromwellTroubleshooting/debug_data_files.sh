#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

echo "=== Checking host file system ==="
echo "Contents of $PROJECT_ROOT:"
ls -la "$PROJECT_ROOT"

echo -e "\n=== Checking data directory ==="
if [ -d "$PROJECT_ROOT/data" ]; then
    echo "Data directory exists. Contents:"
    ls -la "$PROJECT_ROOT/data"
    
    echo -e "\n=== Checking for specific required files ==="
    for file in "participant.csv" "participant_environment.csv" "covar.txt" "ukb_hla_v2.txt" "hla_participant.csv" "withdrawals.csv"; do
        if [ -f "$PROJECT_ROOT/data/$file" ]; then
            echo "✓ Found: $file ($(ls -lh "$PROJECT_ROOT/data/$file" | awk '{print $5}'))"
        else
            echo "✗ Missing: $file"
        fi
    done
else
    echo "ERROR: Data directory does not exist at $PROJECT_ROOT/data"
    echo "Please create it and add your data files."
fi

echo -e "\n=== Testing Docker mount ==="
echo "Testing if Docker can see the files..."

# Test with a simple container
if command -v docker &> /dev/null; then
    docker run --rm \
    -v "$PROJECT_ROOT:/prsInteractive" \
    ubuntu:latest \
    /bin/bash -c "
            echo 'Inside container:';
            echo 'Contents of /prsInteractive:';
            ls -la /prsInteractive/;
            echo '';
            echo 'Contents of /prsInteractive/data:';
            ls -la /prsInteractive/data/;
            echo '';
            echo 'Checking for participant.csv:';
            if [ -f '/prsInteractive/data/participant.csv' ]; then
                echo '✓ participant.csv found in container';
                ls -la /prsInteractive/data/participant.csv;
            else
                echo '✗ participant.csv NOT found in container';
            fi
        "
else
    echo "Docker not available for testing"
fi