#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Quick Cromwell Status Check ==="

# Check if Docker containers are running
echo "1. Docker containers:"
docker ps --filter "ancestor=ukb-base:V1" --format "table {{.Names}}\t{{.Status}}\t{{.Command}}" || echo "No ukb-base containers running"

echo ""
echo "2. Recent Docker container logs:"
CONTAINER_ID=$(docker ps -q --filter "ancestor=ukb-base:V1" | head -1)
if [ -n "$CONTAINER_ID" ]; then
    echo "Container $CONTAINER_ID is running. Last 5 log lines:"
    docker logs --tail 5 "$CONTAINER_ID" 2>&1
else
    echo "No ukb-base containers currently running"
fi

echo ""
echo "3. Cromwell execution progress:"
if [ -d "$PROJECT_ROOT/cromwell-executions" ]; then
    echo "Execution directories:"
    find "$PROJECT_ROOT/cromwell-executions" -maxdepth 3 -type d | head -10
    
    # Look for recent stdout/stderr
    RECENT_STDOUT=$(find "$PROJECT_ROOT/cromwell-executions" -name "stdout" -type f -exec stat -f "%m %N" {} \; 2>/dev/null | sort -nr | head -1 | cut -d' ' -f2-)
    if [ -n "$RECENT_STDOUT" ]; then
        echo ""
        echo "Most recent stdout:"
        tail -5 "$RECENT_STDOUT"
    fi
    
    RECENT_STDERR=$(find "$PROJECT_ROOT/cromwell-executions" -name "stderr" -type f -exec stat -f "%m %N" {} \; 2>/dev/null | sort -nr | head -1 | cut -d' ' -f2-)
    if [ -n "$RECENT_STDERR" ]; then
        echo ""
        echo "Most recent stderr:"
        tail -5 "$RECENT_STDERR"
    fi
else
    echo "No cromwell-executions directory yet"
fi

echo ""
echo "4. System resources:"
echo "CPU usage:"
top -l 1 -n 0 | grep "CPU usage" || echo "CPU info not available"
echo "Memory usage:"
top -l 1 -n 0 | grep "PhysMem" || echo "Memory info not available"

echo ""
echo "=== Status Summary ==="
if [ -n "$CONTAINER_ID" ]; then
    echo "✅ Docker container is running - workflow is active"
    echo "💡 This is normal - the container may take several minutes to complete"
    echo "💡 Check for file access issues in the container logs above"
else
    echo "⏳ No active container - may be between tasks or completed"
fi

echo ""
echo "To monitor continuously: ./monitor_cromwell.sh"
echo "To check specific logs: ls -la cromwell-executions/*/call-*/execution/"