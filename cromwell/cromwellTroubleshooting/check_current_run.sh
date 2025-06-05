#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Checking Current Cromwell Run Status ==="

# Check if Cromwell is still running
echo "1. Cromwell processes:"
ps aux | grep cromwell | grep -v grep || echo "No Cromwell processes running"

echo ""
echo "2. Current Docker containers:"
docker ps --format "table {{.Names}}\t{{.Image}}\t{{.Status}}\t{{.Command}}" | grep ukb-base || echo "No ukb-base containers running"

echo ""
echo "3. Most recent execution results:"
LATEST_EXEC=$(find "$PROJECT_ROOT/cromwell-executions" -name "execution" -type d 2>/dev/null | sort | tail -1)
if [ -n "$LATEST_EXEC" ]; then
    echo "Latest execution: $LATEST_EXEC"
    
    # Check if it completed
    if [ -f "$LATEST_EXEC/rc" ]; then
        RC=$(cat "$LATEST_EXEC/rc")
        echo "Return code: $RC"
        
        if [ "$RC" = "0" ]; then
            echo "✅ Task completed successfully"
            
            # Check outputs
            if [ -f "$LATEST_EXEC/stdout" ]; then
                echo ""
                echo "Recent stdout (last 10 lines):"
                tail -10 "$LATEST_EXEC/stdout"
            fi
        else
            echo "❌ Task failed with return code: $RC"
            
            if [ -f "$LATEST_EXEC/stderr" ]; then
                echo ""
                echo "Error output:"
                cat "$LATEST_EXEC/stderr"
            fi
        fi
    else
        echo "⏳ Task still running (no return code yet)"
    fi
else
    echo "No execution directories found"
fi

echo ""
echo "4. Results check:"
if [ -d "$PROJECT_ROOT/results" ]; then
    echo "Results directory contents:"
    find "$PROJECT_ROOT/results" -type f | head -10
else
    echo "No results directory yet"
fi

echo ""
echo "5. Recent workflow logs:"
if [ -d "$PROJECT_ROOT/logs" ]; then
    LATEST_LOG=$(ls -t "$PROJECT_ROOT/logs"/workflow.*.log | head -1)
    if [ -n "$LATEST_LOG" ]; then
        echo "Latest workflow log: $LATEST_LOG"
        echo "Last 5 lines:"
        tail -5 "$LATEST_LOG"
    fi
fi

echo ""
echo "=== Recommendations ==="
if ps aux | grep cromwell | grep -v grep >/dev/null; then
    echo "✅ Cromwell is still running"
    echo "💡 Wait for current execution to complete, or check logs for progress"
else
    echo "⚠️ Cromwell is not running"
    echo "💡 Either it completed, or it crashed"
    echo "💡 Check the latest logs and execution directories above"
fi

echo ""
echo "To restart with the improved WDL:"
echo "1. Stop current Cromwell (Ctrl+C if running)"
echo "2. Replace your WDL file with the fixed version"
echo "3. Run: ./run_with_docker.sh"

