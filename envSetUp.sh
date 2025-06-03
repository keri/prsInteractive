#!/bin/bash

#set up environment variables for use in workflow
platform=${5:-"local"}  # default to local
pheno=$1
icd10=$2
phenoStr=$3
n=$4 

#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"
#n=40

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "project root is : ${PROJECT_ROOT}"

# Export environment variable for Docker mounting
export PRS_INTERACTIVE_HOME="$PROJECT_ROOT"

# Set base directories
#make a folder inside root data folder for each phenotype

# Create necessary directories
mkdir -p "$PROJECT_ROOT/results/$pheno"
mkdir -p "$PROJECT_ROOT/logs"
mkdir -p "$PROJECT_ROOT/logs/call-logs"
mkdir -p "$PROJECT_ROOT/cromwell-executions"

# Create config directory
mkdir -p "$PROJECT_ROOT/config"

# Create base configuration
cat > "$PROJECT_ROOT/config/default.config" << 'EOF'
# Multi-platform configuration
backend {
    default = "Local"  # Can be overridden at runtime
    
    providers {
        # Local development/testing
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
                    docker run \
                        --rm \
                        --cpus="${cpu}" \
                        --memory="${memory}" \
                        -v ${cwd}:${docker_cwd} \
                        -v ${PRS_INTERACTIVE_HOME}/data:/data:ro \
                        -i ${docker} \
                        /bin/bash < ${script}
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
                    sbatch \
                        --wait \
                        --partition=${partition} \
                        --cpus-per-task=${cpu} \
                        --mem=${memory} \
                        --time=${time_minutes} \
                        --wrap "/bin/bash ${script}"
                """
                
                submit-docker = """
                    sbatch \
                        --wait \
                        --partition=${partition} \
                        --cpus-per-task=${cpu} \
                        --mem=${memory} \
                        --time=${time_minutes} \
                        --wrap "singularity exec --bind ${cwd}:${docker_cwd} docker://${docker} /bin/bash ${script}"
                """
                
                kill = "scancel ${job_id}"
                check-alive = "squeue -j ${job_id}"
                job-id-regex = "Submitted batch job (\\d+)"
            }
        }
        
        # DNAnexus backend
        DNAnexus {
            actor-factory = "cromwell.backend.impl.dnanexus.DxBackendLifecycleActorFactory"
            config {
                default-runtime-attributes {
                    dx_instance_type: "mem1_ssd1_v2_x4"
                    dx_timeout: "1H"
                }
                
                runtime-attributes = """
                    String? dx_instance_type
                    String? dx_timeout
                    Int? dx_restart_max = 1
                """
            }
        }
    }
}

# Optional: Database configuration for call caching across platforms
database {
    profile = "slick.jdbc.HsqldbProfile$"
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
EOF


echo "[WORKFLOW] DATA_PATH is set to: $PROJECT_ROOT/data"
echo "[WORKFLOW] RESULTS_PATH is set to: $PROJECT_ROOT/results"
echo "[WORKFLOW] Scripts directory: $PROJECT_ROOT/scripts"
echo "CONDA ENV BEING ACTIVATED ...: $PROJECT_ROOT/ukb_en"
echo "[WORKFLOW] HPC_DIR IS SET TO: $PROJECT_ROOT/hpc"

bash "$PROJECT_ROOT/update_pipeline_inputs.sh" "$pheno" "$icd10" "$phenoStr" $n



