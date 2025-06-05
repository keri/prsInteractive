#!/bin/bash

# Quick fix for memory format issue in Cromwell configuration
set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

echo "=== Fixing Memory Format Issues ==="

# Fix the Cromwell configuration
echo "Updating Cromwell configuration..."
cat > config/simple.config << 'EOF'
backend {
    default = "Local"
    
    providers {
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                runtime-attributes = """
                    String? docker
                    Int cpu = 1
                    String memory = "2GB"
                    String disk = "10GB"
                    Int? time_minutes = 60
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

# Fix the WDL file
echo "Updating WDL memory specifications..."
if [ -f "workflows/prsInteractivePipeline.DataProcessing.wdl" ]; then
    # Create backup
    cp workflows/prsInteractivePipeline.DataProcessing.wdl workflows/prsInteractivePipeline.DataProcessing.wdl.backup
    
    # Fix memory format (remove spaces, use consistent format)
    sed -i.bak 's/memory: "2 GB"/memory: "2GB"/g' workflows/prsInteractivePipeline.DataProcessing.wdl
    sed -i.bak 's/memory: "50 GB"/memory: "50GB"/g' workflows/prsInteractivePipeline.DataProcessing.wdl
    
    echo "✅ WDL file updated"
    
    # Show what was changed
    echo "Memory specifications updated:"
    grep -n "memory:" workflows/prsInteractivePipeline.DataProcessing.wdl || echo "No memory specifications found"
else
    echo "❌ WDL file not found at workflows/prsInteractivePipeline.DataProcessing.wdl"
fi

echo "✅ Memory format fixes applied!"
echo ""
echo "Now run:"
echo "java -Dconfig.file=config/simple.config -jar cromwell.jar run workflows/prsInteractivePipeline.DataProcessing.wdl --inputs workflows/pipelineInputs.json --options workflows/options.json"