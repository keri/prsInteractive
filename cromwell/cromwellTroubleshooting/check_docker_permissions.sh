#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Docker Permission and Mount Diagnostics ==="
echo "Project root: $PROJECT_ROOT"
echo "Current user: $(whoami)"
echo "Docker version: $(docker --version)"

echo ""
echo "=== Checking file permissions on host ==="
ls -la "$PROJECT_ROOT"
ls -la "$PROJECT_ROOT/data"

echo ""
echo "=== Checking if Docker can access the directory ==="
# Test basic volume mounting
docker run --rm -v "$PROJECT_ROOT:/test" ubuntu:latest ls -la /test

echo ""
echo "=== Testing with absolute paths ==="
# Get the absolute path
ABS_PATH=$(realpath "$PROJECT_ROOT")
echo "Absolute path: $ABS_PATH"

docker run --rm \
-v "$ABS_PATH:/prsInteractive" \
-v "$ABS_PATH/data:/data:ro" \
ubuntu:latest \
/bin/bash -c "
        echo 'Testing with absolute paths:';
        echo 'Contents of /prsInteractive:';
        ls -la /prsInteractive/ | head -10;
        echo '';
        echo 'Contents of /data:';
        ls -la /data/ | head -10;
        echo '';
        echo 'Looking for participant.csv:';
        find /prsInteractive -name 'participant.csv' -type f;
        find /data -name 'participant.csv' -type f;
    "

echo ""
echo "=== Testing with different mount syntax ==="
# Try without :ro flag
docker run --rm \
-v "$ABS_PATH:/prsInteractive" \
-v "$ABS_PATH/data:/data" \
ubuntu:latest \
/bin/bash -c "
        echo 'Testing without readonly flag:';
        ls -la /data/participant.csv 2>/dev/null && echo 'SUCCESS: Can access participant.csv' || echo 'FAILED: Cannot access participant.csv';
    "

echo ""
echo "=== Checking Docker Desktop file sharing (macOS/Windows) ==="
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Running on macOS - check Docker Desktop file sharing settings"
    echo "Make sure $PROJECT_ROOT is in Docker Desktop > Preferences > Resources > File Sharing"
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "Running on Windows - check Docker Desktop file sharing settings"
    echo "Make sure the drive containing $PROJECT_ROOT is shared in Docker Desktop"
else
    echo "Running on Linux - file sharing should work by default"
fi

echo ""
echo "=== SELinux check (Linux only) ==="
if command -v getenforce &> /dev/null; then
    echo "SELinux status: $(getenforce)"
    echo "If SELinux is enforcing, you may need to add :Z flag to volume mounts"
fi

