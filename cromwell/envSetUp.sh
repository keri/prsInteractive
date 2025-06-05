#!/bin/bash

#set up environment variables for use in workflow
platform=${5:-"local"}  # default to local
pheno=$1
icd10=$2
phenoStr=$3
n=$4 

set -e

# Generate a fixed Cromwell configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ABS_PROJECT_ROOT=$(realpath "$PROJECT_ROOT")

echo "Creating fixed Cromwell configuration..."
echo "Project root: $ABS_PROJECT_ROOT"

# Export environment variable for Docker mounting
export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"

# Create necessary directories
mkdir -p "$PROJECT_ROOT/results/$pheno"
mkdir -p "$PROJECT_ROOT/logs"
mkdir -p "$PROJECT_ROOT/logs/call-logs"
mkdir -p "$PROJECT_ROOT/cromwell-executions"

# Create config directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/config"

# Generate the fixed default.config
cat > "$PROJECT_ROOT/config/default.config" << EOF
# Multi-platform configuration with absolute paths
backend {
    default = "Local"
    
    providers {
        # Local development/testing with fixed paths
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                runtime-attributes = """
                    String? docker
                    Int cpu = 1
                    String memory = "2 GB"
                    String disk = "10 GB"
                    Int? time_minutes = 60
                """
                
                submit-docker = """
                    docker run \\
                        --rm \\
                        --cpus="\${cpu}" \\
                        --memory="\${memory}" \\
                        -v \${cwd}:\${docker_cwd} \\
                        -v $ABS_PROJECT_ROOT/data:/prsInteractive/data:ro \\
                        -v $ABS_PROJECT_ROOT/results:/prsInteractive/results \\
                        -v $ABS_PROJECT_ROOT/logs:/prsInteractive/logs \\
                        -v $ABS_PROJECT_ROOT/config:/prsInteractive/config:ro \\
                        -w \${docker_cwd} \\
                        \${docker} \\
                        /bin/bash \${script}
                """
                
                concurrent-job-limit = 4
            }
        }
        
        # Alternative local config with different mount strategy
        LocalAlt {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                runtime-attributes = """
                    String? docker
                    Int cpu = 1
                    String memory = "2 GB"
                    String disk = "10 GB"
                    Int? time_minutes = 60
                """
                
                submit-docker = """
                    docker run \\
                        --rm \\
                        --cpus="\${cpu}" \\
                        --memory="\${memory}" \\
                        -v \${cwd}:\${docker_cwd} \\
                        -v $ABS_PROJECT_ROOT:/host_prsInteractive:ro \\
                        -w \${docker_cwd} \\
                        \${docker} \\
                        /bin/bash \${script}
                """
                
                concurrent-job-limit = 4
            }
        }
        
        # HPC/SLURM backend
        SLURM {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                runtime-attributes = """
                    String? docker
                    String? singularity
                    Int cpu = 1
                    String memory = "2 GB"
                    String disk = "10 GB"
                    Int? time_minutes = 60
                    String? partition = "compute"
                    String? account
                    String? qos
                """
                
                submit = """
                    sbatch \\
                        --wait \\
                        --partition=\${partition} \\
                        --cpus-per-task=\${cpu} \\
                        --mem=\${memory} \\
                        --time=\${time_minutes} \\
                        --wrap "/bin/bash \${script}"
                """
                
                submit-docker = """
                    sbatch \\
                        --wait \\
                        --partition=\${partition} \\
                        --cpus-per-task=\${cpu} \\
                        --mem=\${memory} \\
                        --time=\${time_minutes} \\
                        --wrap "singularity exec \\
                            --bind \${cwd}:\${docker_cwd} \\
                            --bind $ABS_PROJECT_ROOT/data:/prsInteractive/data:ro \\
                            --bind $ABS_PROJECT_ROOT/results:/prsInteractive/results \\
                            docker://\${docker} /bin/bash \${script}"
                """
                
                kill = "scancel \${job_id}"
                check-alive = "squeue -j \${job_id}"
                job-id-regex = "Submitted batch job (\\\\d+)"
            }
        }
    }
}

database {
    profile = "slick.jdbc.HsqldbProfile\$"
    db {
        driver = "org.hsqldb.jdbcDriver"
        url = "jdbc:hsqldb:mem:cromwell-db;shutdown=true;hsqldb.tx=mvcc"
        connectionTimeout = 1000
    }
}
EOF

# Create platform-specific configs
cat > "$PROJECT_ROOT/config/local.config" << 'EOF'
include "default.config"
backend.default = "Local"
EOF

cat > "$PROJECT_ROOT/config/local-alt.config" << 'EOF'
include "default.config"
backend.default = "LocalAlt"
EOF

cat > "$PROJECT_ROOT/config/hpc.config" << 'EOF'
include "default.config"
backend.default = "SLURM"
EOF

cat > "$PROJECT_ROOT/config/dnanexus.config" << 'EOF'
include "default.config"
backend.default = "DNAnexus"
EOF

# Add this after creating the configs
case $platform in
    "local")
        CONFIG_FILE="$PROJECT_ROOT/config/local.config"
        echo "Using local configuration"
    ;;
    "local-alt")
        CONFIG_FILE="$PROJECT_ROOT/config/local-alt.config"
        echo "Using alternative local configuration"
    ;;
    "hpc")
        CONFIG_FILE="$PROJECT_ROOT/config/hpc.config"
        echo "Using HPC/SLURM configuration"
    ;;
    "dnanexus")
        CONFIG_FILE="$PROJECT_ROOT/config/dnanexus.config"
        echo "Using DNAnexus configuration"
    ;;
    *)
        echo "Unknown platform: $platform. Using local."
        CONFIG_FILE="$PROJECT_ROOT/config/local.config"
    ;;
esac

# Export for later use
export CROMWELL_CONFIG="$CONFIG_FILE"

echo "Fixed configuration created!"
echo "Config file: $CONFIG_FILE"

# Create a .env file for future reference
cat > "$PROJECT_ROOT/.env" << EOF
# prsInteractive Pipeline Environment
PRS_INTERACTIVE_HOME=$PROJECT_ROOT
CROMWELL_ROOT=$PROJECT_ROOT/cromwell-executions
RESULTS_PATH=$PROJECT_ROOT/results
LOGS_DIR=$PROJECT_ROOT/logs
WORKFLOW_DIR=$PROJECT_ROOT/workflow
SCRIPTS_DIR=$PROJECT_ROOT/scripts
DATA_PATH=$PROJECT_ROOT/data
HPC_DIR=$PROJECT_ROOT/hpc
ENV_PATH=$PROJECT_ROOT/ukb_env
PLATFORM=$platform
CROMWELL_CONFIG=$CONFIG_FILE
EOF

echo "[WORKFLOW] DATA_PATH is set to: $PROJECT_ROOT/data"
echo "[WORKFLOW] RESULTS_PATH is set to: $PROJECT_ROOT/results"
echo "[WORKFLOW] Scripts directory: $PROJECT_ROOT/scripts"
echo "CONDA ENV BEING ACTIVATED ...: $PROJECT_ROOT/ukb_env"
echo "[WORKFLOW] HPC_DIR IS SET TO: $PROJECT_ROOT/hpc"
echo "[PLATFORM] PLATFORM USED TO RUN ANALYSIS: $platform"

# Make sure data directory exists
if [ ! -d "$PROJECT_ROOT/data" ]; then
    echo "WARNING: Data directory $PROJECT_ROOT/data does not exist!"
    echo "Please create it and place your data files there."
fi

# Check for required data files
echo "Checking for required data files..."
for file in "participant.csv" "participant_environment.csv" "covar.txt" "ukb_hla_v2.txt" "hla_participant.csv" "withdrawals.csv"; do
    if [ ! -f "$PROJECT_ROOT/data/$file" ]; then
        echo "WARNING: Required file $PROJECT_ROOT/data/$file not found!"
    else
        echo "✓ Found: $file"
    fi
done

bash "$PROJECT_ROOT/update_pipeline_inputs.sh" "$pheno" "$icd10" "$phenoStr" $n "$platform"


