#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Current Pipeline Configuration Check ==="
echo "Project root: $PROJECT_ROOT"

echo ""
echo "=== Checking pipelineInputs.json ==="
if [ -f "$PROJECT_ROOT/workflows/pipelineInputs.json" ]; then
    echo "✅ pipelineInputs.json exists"
    echo ""
    echo "Current configuration:"
    cat "$PROJECT_ROOT/workflows/pipelineInputs.json" | python3 -m json.tool 2>/dev/null || cat "$PROJECT_ROOT/workflows/pipelineInputs.json"
    
    echo ""
    echo "=== Key paths being used ==="
    grep -E "(config_file|data|scripts)" "$PROJECT_ROOT/workflows/pipelineInputs.json" || echo "No key paths found"
else
    echo "❌ pipelineInputs.json not found"
    echo "Expected location: $PROJECT_ROOT/workflows/pipelineInputs.json"
    ls -la "$PROJECT_ROOT/workflows/" 2>/dev/null || echo "workflows directory doesn't exist"
fi

echo ""
echo "=== Checking Cromwell config file ==="
if [ -f "$PROJECT_ROOT/.env" ]; then
    source "$PROJECT_ROOT/.env"
    echo "Environment file loaded"
    echo "CROMWELL_CONFIG: $CROMWELL_CONFIG"
    
    if [ -f "$CROMWELL_CONFIG" ]; then
        echo "✅ Cromwell config file exists: $CROMWELL_CONFIG"
    else
        echo "❌ Cromwell config file missing: $CROMWELL_CONFIG"
    fi
else
    echo "❌ .env file not found"
fi

echo ""
echo "=== Checking scripts archive ==="
if [ -f "$PROJECT_ROOT/scripts.tar.gz" ]; then
    echo "✅ scripts.tar.gz exists"
    echo "Contents:"
    tar -tzf "$PROJECT_ROOT/scripts.tar.gz" | head -10
else
    echo "❌ scripts.tar.gz missing"
    echo "Creating scripts archive..."
    cd "$PROJECT_ROOT"
    if [ -d "scripts" ]; then
        tar -czf scripts.tar.gz scripts/
        echo "✅ Created scripts.tar.gz"
    else
        echo "❌ scripts directory doesn't exist"
    fi
fi

echo ""
echo "=== Checking Python script ==="
if [ -d "$PROJECT_ROOT/scripts" ]; then
    find "$PROJECT_ROOT/scripts" -name "create_pheno_train_test_split.py" -type f
    if [ $? -eq 0 ]; then
        echo "✅ Python script found"
    else
        echo "❌ create_pheno_train_test_split.py not found in scripts directory"
        echo "Scripts directory contents:"
        ls -la "$PROJECT_ROOT/scripts/"
    fi
else
    echo "❌ scripts directory doesn't exist"
fi

echo ""
echo "=== Most recent Cromwell error log ==="
if [ -d "$PROJECT_ROOT/cromwell-executions" ]; then
    # Find the most recent cromwell execution
    LATEST_EXEC=$(find "$PROJECT_ROOT/cromwell-executions" -name "stderr" -type f -exec stat -f "%m %N" {} \; 2>/dev/null | sort -nr | head -1 | cut -d' ' -f2-)
    if [ -n "$LATEST_EXEC" ]; then
        echo "Latest error log: $LATEST_EXEC"
        echo "Last 20 lines:"
        tail -20 "$LATEST_EXEC"
    else
        echo "No recent error logs found"
    fi
else
    echo "No cromwell-executions directory found"
fi

