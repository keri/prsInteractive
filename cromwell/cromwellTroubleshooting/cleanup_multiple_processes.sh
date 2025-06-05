#!/bin/bash

echo "=== Cleaning Up Multiple Cromwell Processes ==="

echo "1. Current Cromwell processes:"
ps aux | grep cromwell | grep -v grep

echo ""
echo "2. Stopping all Cromwell processes..."

# Kill all java processes running cromwell
pkill -f "cromwell.*jar"

echo "Waiting 5 seconds for processes to stop..."
sleep 5

echo ""
echo "3. Checking if any processes remain:"
REMAINING=$(ps aux | grep cromwell | grep -v grep || echo "No cromwell processes found")
if [ "$REMAINING" != "No cromwell processes found" ]; then
    echo "Some processes still running:"
    echo "$REMAINING"
    echo ""
    echo "Force killing remaining processes..."
    pkill -9 -f "cromwell.*jar"
    sleep 2
fi

echo ""
echo "4. Final check:"
ps aux | grep cromwell | grep -v grep || echo "✅ All Cromwell processes stopped"

echo ""
echo "5. Cleaning up any orphaned Docker containers:"
# Check for any ukb-base containers that might be orphaned
ORPHANED=$(docker ps -a --filter "ancestor=ukb-base:V1" --format "{{.ID}}" || echo "")
if [ -n "$ORPHANED" ]; then
    echo "Found orphaned containers, cleaning up:"
    docker rm -f $ORPHANED
else
    echo "✅ No orphaned containers found"
fi

echo ""
echo "6. Cleaning up execution state (optional):"
echo "Do you want to clean up previous execution directories? [y/N]"
read -r response
if [[ "$response" =~ ^[Yy]$ ]]; then
    echo "Cleaning up cromwell-executions..."
    rm -rf cromwell-executions/*/
    echo "✅ Execution directories cleaned"
else
    echo "Keeping existing execution directories"
fi

echo ""
echo "=== Cleanup Complete ==="
echo ""
echo "You can now run a single Cromwell instance:"
echo "  ./run_with_docker.sh"
echo ""
echo "Tips to avoid this in the future:"
echo "1. Always stop previous Cromwell runs (Ctrl+C) before starting new ones"
echo "2. Check 'ps aux | grep cromwell' before starting"
echo "3. Use 'pkill -f cromwell' to stop all instances if needed"

