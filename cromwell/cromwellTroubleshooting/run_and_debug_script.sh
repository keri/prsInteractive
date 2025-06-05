#!/bin/bash

# Complete debug and run script for prsInteractive pipeline
set -e

echo "=== prsInteractive Pipeline Debug and Run Script ==="

# Configuration
PHENO=${1:-"type2Diabetes"}
ICD10=${2:-"E11"}
PHENO_STR=${3:-"type 2 diabetes"}
N=${4:-40}
PLATFORM=${5:-"local"}

# Get project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Project root: $PROJECT_ROOT"

# Step 1: Check Docker installation and permissions
echo ""
echo "=== Step 1: Docker Check ==="
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed or not in PATH"
    exit 1
fi

echo "✅ Docker found: $(docker --version)"

# Test docker permissions
if ! docker ps &> /dev/null; then
    echo "❌ Docker permission denied. Try running with sudo or add user to docker group"
    exit 1
fi

echo "✅ Docker permissions OK"

# Step 2: Check required files
echo ""
echo "=== Step 2: File Structure Check ==="

# Check data files
required_data_files=("participant.csv" "participant_environment.csv" "covar.txt" "hla_participant.csv" "withdrawals.csv" "ukb_hla_v2.txt")
for file in "${required_data_files[@]}"; do
    if [ -f "$PROJECT_ROOT/data/$file" ]; then
        echo "✅ $file found"
    else
        echo "❌ $file missing from $PROJECT_ROOT/data/"
    fi
done

# Check scripts
if [ -f "$PROJECT_ROOT/scripts/create_pheno_train_test_split.py" ]; then
    echo "✅ Python script found"
else
    echo "❌ create_pheno_train_test_split.py missing from scripts/"
fi

# Check if scripts archive exists, if not create it
if [ ! -f "$PROJECT_ROOT/scripts.tar.gz" ]; then
    echo "Creating scripts archive..."
    cd "$PROJECT_ROOT"
    tar -czf scripts.tar.gz -C scripts .
    echo "✅ scripts.tar.gz created"
else
    echo "✅ scripts.tar.gz found"
fi

# Step 3: Build Docker image
echo ""
echo "=== Step 3: Docker Image Build ==="

# Check if Dockerfile exists
if [ ! -f "$PROJECT_ROOT/docker/Dockerfile.base" ]; then
    echo "❌ Dockerfile.base not found in docker/ directory"
    exit 1
fi

# Check if environment.yml exists
if [ ! -f "$PROJECT_ROOT/environment.yml" ]; then
    echo "❌ environment.yml not found"
    exit 1
fi

echo "Building Docker image ukb-base:V1..."
cd "$PROJECT_ROOT"

# Copy necessary files to docker directory
cp environment.yml docker/
mkdir -p docker/config
cp -r config/* docker/config/ 2>/dev/null || echo "No config files to copy"

# Build image
docker build -f docker/Dockerfile.base -t ukb-base:V1 docker/

if [ $? -eq 0 ]; then
    echo "✅ Docker image built successfully"
else
    echo "❌ Docker image build failed"
    exit 1
fi

# Step 4: Test Docker image
echo ""
echo "=== Step 4: Docker Image Test ==="

echo "Testing Docker image..."
docker run --rm ukb-base:V1 python --version
docker run --rm ukb-base:V1 conda list | grep pandas

if [ $? -eq 0 ]; then
    echo "✅ Docker image working correctly"
else
    echo "❌ Docker image test failed"
    exit 1
fi

# Step 5: Create necessary directories
echo ""
echo "=== Step 5: Directory Setup ==="

mkdir -p "$PROJECT_ROOT/results/$PHENO"
mkdir -p "$PROJECT_ROOT/logs"
mkdir -p "$PROJECT_ROOT/cromwell-executions"
mkdir -p "$PROJECT_ROOT/workflows"

echo "✅ Directories created"

# Step 6: Create simplified Cromwell configuration
echo ""
echo "=== Step 6: Cromwell Configuration ==="

cat > "$PROJECT_ROOT/config/simple.config" << 'EOF'
backend {
    default = "Local"
    
    providers {
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                runtime-attributes = """
                    String? docker
                    Int cpu = 1
                    String memory = "2 GB"
                """
                
                submit-docker = """
                    docker run \
                        --rm \
                        --cpus="${cpu}" \
                        --memory="${memory}" \
                        -v ${cwd}:${docker_cwd} \
                        -w ${docker_cwd} \
                        ${docker} \
                        /bin/bash ${script}
                """
                
                concurrent-job-limit = 1
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

echo "✅ Simplified Cromwell config created"

# Step 7: Create pipeline inputs
echo ""
echo "=== Step 7: Pipeline Inputs ==="

cat > "$PROJECT_ROOT/workflows/pipelineInputs.json" << EOF
{
    "prsInteractivePipeline.pheno": "$PHENO",
    "prsInteractivePipeline.icd10": "$ICD10",
    "prsInteractivePipeline.pheno_str": "$PHENO_STR",
    "prsInteractivePipeline.n": $N,
    "prsInteractivePipeline.scripts_archive": "$PROJECT_ROOT/scripts.tar.gz",
    "prsInteractivePipeline.root_directory": "$PROJECT_ROOT",
    "prsInteractivePipeline.platform": "$PLATFORM",
    "prsInteractivePipeline.raw_variant_calls": [
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c1_b0_v2"
    ],
    "prsInteractivePipeline.participant_data": "$PROJECT_ROOT/data/participant.csv",
    "prsInteractivePipeline.participant_environment": "$PROJECT_ROOT/data/participant_environment.csv",
    "prsInteractivePipeline.covar": "$PROJECT_ROOT/data/covar.txt",
    "prsInteractivePipeline.hla_headers": "$PROJECT_ROOT/data/ukb_hla_v2.txt",
    "prsInteractivePipeline.participant_hla": "$PROJECT_ROOT/data/hla_participant.csv",
    "prsInteractivePipeline.withdrawals": "$PROJECT_ROOT/data/withdrawals.csv"
}
EOF

echo "✅ Pipeline inputs created"

# Step 8: Create options file
echo ""
echo "=== Step 8: Options Configuration ==="

cat > "$PROJECT_ROOT/workflows/options.json" << 'EOF'
{
  "write_to_cache": false,
  "read_from_cache": false,
  "final_workflow_outputs_dir": "./results",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "./logs",
  "final_call_logs_dir": "./logs/call-logs",
  "delete_intermediate_output_files": false,
  "workflow_failure_mode": "ContinueWhilePossible",
  "default_runtime_attributes": {
    "docker": "ukb-base:V1"
  },
  "workflow_log_level": "INFO"
}
EOF

echo "✅ Options file created"

# Step 9: Test Python script standalone
echo ""
echo "=== Step 9: Python Script Test ==="

echo "Testing Python script in Docker container..."

# Create a test container to verify the script works
docker run --rm \
    -v "$PROJECT_ROOT/data:/data:ro" \
    -v "$PROJECT_ROOT/results:/results" \
    -v "$PROJECT_ROOT/scripts:/scripts:ro" \
    ukb-base:V1 \
    python /scripts/create_pheno_train_test_split.py \
        --data_path /data \
        --pheno_path /results/test_$PHENO \
        --pheno $PHENO \
        --pheno_str "$PHENO_STR" \
        --icd10 $ICD10

if [ $? -eq 0 ]; then
    echo "✅ Python script test successful"
    echo "Test output files:"
    ls -la "$PROJECT_ROOT/results/test_$PHENO/" || echo "No output directory created"
else
    echo "❌ Python script test failed"
    echo "Check the error messages above for details"
fi

# Step 10: Check Cromwell installation
echo ""
echo "=== Step 10: Cromwell Check ==="

# Check if Cromwell JAR exists
if [ ! -f "$PROJECT_ROOT/cromwell.jar" ]; then
    echo "Cromwell JAR not found. Downloading..."
    cd "$PROJECT_ROOT"
    wget -O cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
    if [ $? -eq 0 ]; then
        echo "✅ Cromwell downloaded successfully"
    else
        echo "❌ Failed to download Cromwell"
        exit 1
    fi
else
    echo "✅ Cromwell JAR found"
fi

# Test Java installation
if ! command -v java &> /dev/null; then
    echo "❌ Java is not installed"
    exit 1
fi

echo "✅ Java found: $(java -version 2>&1 | head -1)"

# Step 11: Run the workflow
echo ""
echo "=== Step 11: Running Workflow ==="

cd "$PROJECT_ROOT"

echo "Starting Cromwell workflow..."
echo "Phenotype: $PHENO"
echo "ICD10: $ICD10"
echo "Phenotype string: $PHENO_STR"

# Clean up previous executions
rm -rf cromwell-executions/*

# Run Cromwell (version 85+ syntax)
java -Dconfig.file=config/simple.config -jar cromwell.jar run \
    workflows/prsInteractivePipeline.DataProcessing.wdl \
    --inputs workflows/pipelineInputs.json \
    --options workflows/options.json \
    2>&1 | tee logs/workflow_$(date +%Y%m%d_%H%M%S).log

CROMWELL_EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=== Step 12: Results Check ==="

if [ $CROMWELL_EXIT_CODE -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    # Check output files
    echo "Output files in results/$PHENO/:"
    ls -la "results/$PHENO/" || echo "No output directory found"
    
    # Check for expected files
    for file in "trainingID.txt" "testID.txt" "holdoutID.txt" "pheno.txt"; do
        if [ -f "results/$PHENO/$file" ]; then
            echo "✅ $file created ($(wc -l < "results/$PHENO/$file") lines)"
        else
            echo "❌ $file missing"
        fi
    done
    
else
    echo "❌ Workflow failed with exit code $CROMWELL_EXIT_CODE"
    echo ""
    echo "Check the log files for details:"
    echo "- Main log: logs/workflow_*.log"
    echo "- Cromwell executions: cromwell-executions/"
    
    # Show recent error logs
    echo ""
    echo "Recent error logs:"
    find cromwell-executions -name "stderr" -exec echo "=== {} ===" \; -exec cat {} \; 2>/dev/null | tail -50
fi

echo ""
echo "=== Debugging Information ==="
echo "Docker images:"
docker images | grep ukb-base

echo ""
echo "Project structure:"
find "$PROJECT_ROOT" -maxdepth 2 -type d | sort

echo ""
echo "Data files:"
ls -la "$PROJECT_ROOT/data/" 2>/dev/null || echo "Data directory not accessible"

echo ""
echo "Results:"
ls -la "$PROJECT_ROOT/results/" 2>/dev/null || echo "Results directory not accessible"

echo ""
echo "=== End of Debug Script ==="