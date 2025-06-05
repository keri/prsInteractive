#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Debugging Cromwell Execution Status ==="

# Check all recent Docker containers (including stopped ones)
echo "1. All recent Docker containers (including stopped):"
docker ps -a --filter "ancestor=ukb-base:V1" --format "table {{.Names}}\t{{.Status}}\t{{.CreatedAt}}\t{{.Command}}" | head -10

echo ""
echo "2. Cromwell execution directories:"
if [ -d "$PROJECT_ROOT/cromwell-executions" ]; then
    echo "Finding execution directories..."
    
    # Find all execution directories
    EXEC_DIRS=$(find "$PROJECT_ROOT/cromwell-executions" -name "execution" -type d 2>/dev/null | sort)
    
    if [ -n "$EXEC_DIRS" ]; then
        echo "Found execution directories:"
        echo "$EXEC_DIRS"
        
        # Check the most recent one
        LATEST_EXEC=$(echo "$EXEC_DIRS" | tail -1)
        echo ""
        echo "=== Latest execution directory: $LATEST_EXEC ==="
        
        if [ -d "$LATEST_EXEC" ]; then
            echo "Contents of execution directory:"
            ls -la "$LATEST_EXEC"
            
            echo ""
            echo "=== Checking execution files ==="
            
            # Check script file
            if [ -f "$LATEST_EXEC/script" ]; then
                echo "Generated script exists. Content:"
                cat "$LATEST_EXEC/script"
            else
                echo "❌ No script file found"
            fi
            
            echo ""
            # Check stdout
            if [ -f "$LATEST_EXEC/stdout" ]; then
                echo "STDOUT content:"
                cat "$LATEST_EXEC/stdout"
            else
                echo "❌ No stdout file found"
            fi
            
            echo ""
            # Check stderr
            if [ -f "$LATEST_EXEC/stderr" ]; then
                echo "STDERR content:"
                cat "$LATEST_EXEC/stderr"
            else
                echo "❌ No stderr file found"
            fi
            
            echo ""
            # Check return code
            if [ -f "$LATEST_EXEC/rc" ]; then
                echo "Return code:"
                cat "$LATEST_EXEC/rc"
            else
                echo "❌ No return code file found (task may still be running)"
            fi
            
        else
            echo "❌ Execution directory not accessible"
        fi
    else
        echo "❌ No execution directories found"
    fi
else
    echo "❌ cromwell-executions directory doesn't exist"
fi

echo ""
echo "3. Checking workflow metadata:"
# Look for workflow logs
if [ -d "$PROJECT_ROOT/logs" ]; then
    echo "Logs directory contents:"
    ls -la "$PROJECT_ROOT/logs"
    
    # Check call logs
    if [ -d "$PROJECT_ROOT/logs/call-logs" ]; then
        echo ""
        echo "Call logs:"
        find "$PROJECT_ROOT/logs/call-logs" -type f | head -5
    fi
fi

echo ""
echo "4. Checking Docker daemon and image:"
echo "Docker daemon status:"
docker info >/dev/null 2>&1 && echo "✅ Docker daemon running" || echo "❌ Docker daemon not running"

echo ""
echo "ukb-base:V1 image status:"
docker images ukb-base:V1 --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

echo ""
echo "5. Recent Docker events:"
echo "Docker events in the last 5 minutes:"
docker events --since="5m" --until="now" | grep -E "(ukb-base|create|start|die)" | tail -10 || echo "No recent ukb-base events"

echo ""
echo "=== Diagnosis ==="
if [ -f "$PROJECT_ROOT/cromwell-executions"/*/*/execution/rc ]; then
    echo "✅ Task appears to have completed (return code file exists)"
elif [ -f "$PROJECT_ROOT/cromwell-executions"/*/*/execution/script ]; then
    echo "⏳ Task script generated but no return code yet - may still be running"
else
    echo "❌ No execution artifacts found - task may not have started properly"
fi

