#!/bin/bash

echo "=== Checking Docker Desktop File Sharing Status ==="

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

# Test if Docker can access the directory without explicit volume mount
echo ""
echo "=== Testing basic Docker access ==="
docker run --rm -v "$PROJECT_ROOT:/test" ubuntu:latest ls -la /test > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Docker can access your project directory"
else
    echo "❌ Docker cannot access your project directory"
fi

# Test specific common paths that should be shared by default
echo ""
echo "=== Testing default Docker Desktop shared paths ==="

for path in "/Users" "/tmp" "/private" "/var/folders"; do
    if [ -d "$path" ]; then
        docker run --rm -v "$path:/test" ubuntu:latest ls /test > /dev/null 2>&1
        if [ $? -eq 0 ]; then
            echo "✅ $path is accessible"
        else
            echo "❌ $path is not accessible"
        fi
    fi
done

# Check if we can access the parent directories
echo ""
echo "=== Testing parent directory access ==="
docker run --rm -v "/Users/kerimulterer:/test" ubuntu:latest ls -la /test > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ /Users/kerimulterer is accessible"
    # Check if our project is visible
    docker run --rm -v "/Users/kerimulterer:/test" ubuntu:latest ls -la /test/prsInteractive > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "✅ prsInteractive directory is visible through parent mount"
    else
        echo "❌ prsInteractive directory is not visible through parent mount"
    fi
else
    echo "❌ /Users/kerimulterer is not accessible"
fi

echo ""
echo "=== Docker Desktop Information ==="
echo "Docker version: $(docker --version)"
echo "Docker context: $(docker context show 2>/dev/null || echo 'default')"

# Check Docker Desktop VM details
echo ""
echo "=== Checking Docker Desktop VM ==="
docker run --rm --privileged ubuntu:latest /bin/bash -c "
    echo 'Inside Docker VM:';
    echo 'Mount points:';
    mount | grep -E '(Users|home)' || echo 'No Users/home mounts found';
    echo '';
    echo 'Available directories under /':
    ls -la / | grep -E '(Users|home|mnt|media)' || echo 'No relevant mount points found';
"