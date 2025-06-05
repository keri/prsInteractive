#!/bin/bash

# Test script to verify Docker mounting works correctly
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Testing Docker volume mounting..."
echo "Project root: $PROJECT_ROOT"

# Test the Docker command that Cromwell would use
docker run \
    --rm \
    -v $PROJECT_ROOT/data:/data:ro \
    -v $PROJECT_ROOT:/prsInteractive \
    ubuntu:latest \
    /bin/bash -c "
        echo 'Inside container:';
        echo 'Current working directory: \$(pwd)';
        echo 'Checking /data mount:';
        ls -la /data/ || echo '/data mount failed';
        echo 'Checking /prsInteractive mount:';
        ls -la /prsInteractive/ || echo '/prsInteractive mount failed';
        echo 'Checking for specific files:';
        ls -la /prsInteractive/data/participant.csv || echo 'participant.csv not found';
        ls -la /prsInteractive/data/participant_environment.csv || echo 'participant_environment.csv not found';
    "

echo "Test completed."