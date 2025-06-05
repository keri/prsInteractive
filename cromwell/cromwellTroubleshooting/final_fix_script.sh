#!/bin/bash

# Final fix for Cromwell configuration issues
set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

echo "=== Final Cromwell Configuration Fix ==="

# Kill any existing Cromwell processes
pkill -f cromwell.jar || true
sleep 2

# Check if the ultra-simple config exists
if [ ! -f "config/ultra-simple.config" ]; then
    echo "Creating ultra-simple config..."
    mkdir -p config
    cat > config/ultra-simple.config << 'EOF'
backend {
    default = "Local"
    
    providers {
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                submit-docker = """
                    docker run \
                        --rm \
                        -v ${cwd}:${docker_cwd} \
                        -w ${docker_cwd} \
                        ${docker} \
                        /bin/bash ${script}
                """
            }
        }
    }
}

database {
    profile = "slick.jdbc.HsqldbProfile$"
    db {
        driver = "org.hsqldb.jdbcDriver"
        url = "jdbc:hsqldb:mem:cromwell-db;shutdown=true;hsqldb.tx=mvcc"
        connectionTimeout = 1000
    }
}
EOF
fi


echo "Trying Cromwell with ultra-simple configuration..."
echo "Command: java -Dconfig.file=config/ultra-simple.config -jar cromwell.jar run workflows/prsInteractivePipeline.DataProcessing.wdl --inputs workflows/pipelineInputs.json --options workflows/options.json"
echo ""
echo "Starting Cromwell... (Press Ctrl+C if it hangs with memory errors)"

# Start Cromwell in background and capture PID
java -Dconfig.file=config/ultra-simple.config -jar cromwell.jar run \
workflows/prsInteractivePipeline.DataProcessing.wdl \
--inputs workflows/pipelineInputs.json \
--options workflows/options.json &

CROMWELL_PID=$!

# Wait a bit and check if it's still running
sleep 10

if kill -0 $CROMWELL_PID 2>/dev/null; then
    echo "Cromwell is running - waiting for completion..."
    echo "If you see memory errors above, press Ctrl+C and we'll use direct Docker execution"
    wait $CROMWELL_PID
    CROMWELL_EXIT_CODE=$?
    
    if [ $CROMWELL_EXIT_CODE -eq 0 ]; then
        echo "✅ Cromwell completed successfully!"
        echo "Checking results..."
        ls -la results/type2Diabetes/ 2>/dev/null || echo "No results directory found"
        exit 0
    else
        echo "❌ Cromwell failed with exit code $CROMWELL_EXIT_CODE"
    fi
else
    echo "❌ Cromwell process died quickly"
fi

echo ""
echo "=== Option 2: Direct Docker Execution (Bypass Cromwell) ==="

echo "Running pipeline directly with Docker..."

# Ensure directories exist
mkdir -p results/type2Diabetes

# Parameters
PHENO="type2Diabetes"
ICD10="E11"
PHENO_STR="type 2 diabetes"

echo "Parameters:"
echo "  Phenotype: $PHENO"
echo "  ICD10: $ICD10"
echo "  Phenotype string: $PHENO_STR"
echo ""

# Check if required files exist
echo "Checking required files..."
if [ ! -f "data/participant.csv" ]; then
    echo "❌ data/participant.csv not found"
    exit 1
fi
echo "✅ participant.csv found"

if [ ! -f "data/participant_environment.csv" ]; then
    echo "❌ data/participant_environment.csv not found"
    exit 1
fi
echo "✅ participant_environment.csv found"

if [ ! -f "scripts/create_pheno_train_test_split.py" ]; then
    echo "❌ scripts/create_pheno_train_test_split.py not found"
    exit 1
fi
echo "✅ Python script found"

# Check Docker image
echo "Checking Docker image..."
if ! docker images | grep -q ukb-base; then
    echo "❌ ukb-base:V1 Docker image not found. Please build it first."
    exit 1
fi
echo "✅ Docker image found"

echo ""
echo "Running Python script in Docker container..."

docker run --rm \
-v "$(pwd)/data:/data:ro" \
-v "$(pwd)/results:/results" \
-v "$(pwd)/scripts:/scripts:ro" \
ukb-base:V1 \
bash -c "
        echo 'Container started successfully'
        echo 'Working directory:' \$(pwd)
        echo 'Available Python:' \$(which python)
        echo 'Python version:' \$(python --version)
        echo 'Conda environment:' \$CONDA_DEFAULT_ENV
        
        echo 'Files in /data:'
        ls -la /data/
        
        echo 'Files in /scripts:'
        ls -la /scripts/
        
        echo 'Running phenotype creation script...'
        cd /scripts
        python create_pheno_train_test_split.py \
            --data_path /data \
            --pheno_path /results/$PHENO \
            --pheno $PHENO \
            --pheno_str '$PHENO_STR' \
            --icd10 $ICD10
    "

DOCKER_EXIT_CODE=$?

if [ $DOCKER_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✅ Direct Docker execution completed successfully!"
    echo ""
    echo "=== Results ==="
    if [ -d "results/$PHENO" ]; then
        echo "Output files in results/$PHENO/:"
        ls -la "results/$PHENO/"
        
        echo ""
        echo "File contents summary:"
        for file in "results/$PHENO"/*.txt; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                lines=$(wc -l < "$file" 2>/dev/null || echo "0")
                echo "  $filename: $lines lines"
                if [ "$lines" -gt 0 ] && [ "$lines" -lt 10 ]; then
                    echo "    First few lines:"
                    head -3 "$file" 2>/dev/null | sed 's/^/      /'
                fi
            fi
        done
    else
        echo "❌ No results directory created"
    fi
else
    echo "❌ Direct Docker execution failed with exit code $DOCKER_EXIT_CODE"
    exit 1
fi

echo ""
echo "=== Pipeline Summary ==="
echo "✅ Pipeline completed using direct Docker execution"
echo "📁 Results location: results/$PHENO/"
echo "📊 Phenotype: $PHENO (ICD10: $ICD10)"
echo ""
echo "Expected output files:"
echo "  - trainingID.txt (training set participant IDs)"
echo "  - testID.txt (test set participant IDs)"  
echo "  - holdoutID.txt (holdout set participant IDs)"
echo "  - pheno.txt (all participants with phenotype labels)"