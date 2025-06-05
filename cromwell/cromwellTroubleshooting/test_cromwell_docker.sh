#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

# Test the exact Docker command that Cromwell would use
echo "=== Testing Docker mount command ==="

# Create a simple test script
cat > test_script.sh << 'EOF'
echo "=== Inside Docker Container ==="
echo "Current working directory: $(pwd)"
echo "Current user: $(whoami)"
echo "Environment variables:"
env | grep -E "(HOME|USER|PWD)" || true

echo ""
echo "=== Checking mount points ==="
echo "Root directory contents:"
ls -la / | head -20

echo ""
echo "=== Checking /prsInteractive mount ==="
if [ -d "/prsInteractive" ]; then
    echo "/prsInteractive exists"
    ls -la /prsInteractive/
    
    echo ""
    echo "=== Checking /prsInteractive/data ==="
    if [ -d "/prsInteractive/data" ]; then
        echo "/prsInteractive/data exists"
        ls -la /prsInteractive/data/
        
        echo ""
        echo "=== Checking specific files ==="
        for file in participant.csv participant_environment.csv; do
            if [ -f "/prsInteractive/data/$file" ]; then
                echo "✓ Found: $file"
                ls -la "/prsInteractive/data/$file"
            else
                echo "✗ Missing: $file"
            fi
        done
    else
        echo "✗ /prsInteractive/data does not exist"
    fi
else
    echo "✗ /prsInteractive mount point does not exist"
fi

echo ""
echo "=== Checking /data mount ==="
if [ -d "/data" ]; then
    echo "/data exists"
    ls -la /data/
else
    echo "✗ /data mount point does not exist"
fi
EOF

chmod +x test_script.sh

# Test with the ukb-base:V1 image (if available)
echo "Testing with ukb-base:V1 image..."
if docker image inspect ukb-base:V1 >/dev/null 2>&1; then
    docker run \
    --rm \
    -v "$PROJECT_ROOT:/prsInteractive" \
    -v "$PROJECT_ROOT/data:/data:ro" \
    -v "$(pwd)/test_script.sh:/test_script.sh" \
    ukb-base:V1 \
    /bin/bash /test_script.sh
else
    echo "ukb-base:V1 image not found. Testing with ubuntu..."
    docker run \
    --rm \
    -v "$PROJECT_ROOT:/prsInteractive" \
    -v "$PROJECT_ROOT/data:/data:ro" \
    -v "$(pwd)/test_script.sh:/test_script.sh" \
    ubuntu:latest \
    /bin/bash /test_script.sh
fi

# Clean up
rm -f test_script.sh

echo ""
echo "=== Host verification ==="
echo "Confirming files exist on host:"
ls -la "$PROJECT_ROOT/data/participant.csv"
ls -la "$PROJECT_ROOT/data/participant_environment.csv"

