#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ABS_PROJECT_ROOT=$(realpath "$PROJECT_ROOT")

echo "=== Creating Docker-Forced Configuration ==="

# Create a configuration that forces Docker usage
cat > "$PROJECT_ROOT/config/docker-forced.config" << EOF
# Force Docker configuration for Cromwell 90
include required(classpath("application"))

backend {
    default = "LocalDocker"
    
    providers {
        LocalDocker {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                run-in-background = true
                
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
                        --cpus=\${cpu} \\
                        --memory=\${memory} \\
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

echo "✅ Created docker-forced config: $PROJECT_ROOT/config/docker-forced.config"

# Also create a simpler version without the include
cat > "$PROJECT_ROOT/config/simple-docker.config" << EOF
# Simple Docker configuration for Cromwell 90
backend {
    default = "LocalDocker"
    
    providers {
        LocalDocker {
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
                        --cpus=${cpu} \
                        --memory=${memory}m \
                        -v ${cwd}:\${docker_cwd} \
                        -v $ABS_PROJECT_ROOT/data:/prsInteractive/data:ro \
                        -v $ABS_PROJECT_ROOT/results:/prsInteractive/results \
                        -w ${docker_cwd} \
                        ${docker} \
                        /bin/bash ${script}
                """
                
                concurrent-job-limit = 4
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

echo "✅ Created simple-docker config: $PROJECT_ROOT/config/simple-docker.config"

# Test the configuration
echo ""
echo "=== Testing Docker Command Manually ==="
docker run --rm \
-v "$ABS_PROJECT_ROOT/data:/prsInteractive/data:ro" \
-v "$ABS_PROJECT_ROOT/results:/prsInteractive/results" \
ukb-base:V1 \
/bin/bash -c "ls -la /prsInteractive/data/participant.csv && echo 'SUCCESS: Docker mount working'"

echo ""
echo "=== Creating Updated Run Script ==="
cat > "$PROJECT_ROOT/run_with_docker.sh" << RUNEOF
#!/bin/bash

PROJECT_ROOT="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
cd "\$PROJECT_ROOT"

# Load environment
source .env

CROMWELL_JAR="\${CROMWELL_JAR:-\$HOME/cromwell/cromwell-90.jar}"

echo "=== Running with Docker-Forced Configuration ==="
echo "Config: \$PROJECT_ROOT/config/simple-docker.config"
echo "Command:"
echo "java -Dconfig.file=\$PROJECT_ROOT/config/simple-docker.config -jar \$CROMWELL_JAR run workflows/prsInteractivePipeline.DataProcessing.wdl -i workflows/pipelineInputs.json -o workflows/options.json"
echo ""

java -Dconfig.file="\$PROJECT_ROOT/config/simple-docker.config" -jar "\$CROMWELL_JAR" run \\
    workflows/prsInteractivePipeline.wdl \\
    -i workflows/pipelineInputs.json \\
    -o workflows/options.json
RUNEOF

chmod +x "$PROJECT_ROOT/run_with_docker.sh"

echo "✅ Created run_with_docker.sh"

echo ""
echo "=== What Changed ==="
echo "1. ❌ Old config used 'Local' backend (no Docker)"
echo "2. ✅ New config uses 'LocalDocker' backend (forces Docker)"
echo "3. ✅ Simplified configuration without problematic includes"
echo "4. ✅ Direct absolute paths to your data"

echo ""
echo "=== Next Steps ==="
echo "1. Test with: ./run_with_docker.sh"
echo "2. If that works, update your main script to use the simple-docker config"
echo ""
echo "This should force Cromwell to use Docker containers with proper volume mounts!"

