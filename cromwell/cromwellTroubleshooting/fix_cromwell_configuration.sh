#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Fixing Cromwell Configuration Issue ==="

# 1. First, let's check what config is actually being used
echo "Current pipelineInputs.json config_file:"
grep "config_file" "$PROJECT_ROOT/workflows/pipelineInputs.json"

CONFIG_PATH=$(grep "config_file" "$PROJECT_ROOT/workflows/pipelineInputs.json" | cut -d'"' -f4)
echo "Config path from pipelineInputs.json: $CONFIG_PATH"

if [ -f "$CONFIG_PATH" ]; then
    echo "✅ Config file exists"
    echo "Contents:"
    cat "$CONFIG_PATH"
else
    echo "❌ Config file doesn't exist at: $CONFIG_PATH"
fi

echo ""
echo "=== Regenerating environment setup with proper config ==="

# Run envSetUp to regenerate everything properly
echo "Running: ./envSetUp.sh type2Diabetes E11 \"type 2 diabetes\" 40 local"
./envSetUp.sh type2Diabetes E11 "type 2 diabetes" 40 local

echo ""
echo "=== Verifying the fix ==="

# Check if .env now has CROMWELL_CONFIG
if [ -f "$PROJECT_ROOT/.env" ]; then
    echo "Updated .env file:"
    cat "$PROJECT_ROOT/.env"
    
    source "$PROJECT_ROOT/.env"
    if [ -n "$CROMWELL_CONFIG" ]; then
        echo "✅ CROMWELL_CONFIG is now set: $CROMWELL_CONFIG"
        if [ -f "$CROMWELL_CONFIG" ]; then
            echo "✅ Config file exists"
        else
            echo "❌ Config file still doesn't exist"
        fi
    else
        echo "❌ CROMWELL_CONFIG still empty"
    fi
else
    echo "❌ .env file still doesn't exist"
fi

echo ""
echo "=== Updated pipelineInputs.json ==="
if [ -f "$PROJECT_ROOT/workflows/pipelineInputs.json" ]; then
    echo "New config_file path:"
    grep "config_file" "$PROJECT_ROOT/workflows/pipelineInputs.json"
else
    echo "❌ pipelineInputs.json missing"
fi

