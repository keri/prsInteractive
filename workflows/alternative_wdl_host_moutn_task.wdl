task runTrainTestHoldoutSplitPhenotype {
    input {
        String pheno
        String icd10
        String pheno_str
        String root_directory
        File scripts_archive
    }
    
    command <<<
        set -e
        
        # Extract scripts archive
        echo "Extracting scripts..."
        tar -xzf ~{scripts_archive} || unzip ~{scripts_archive}
        
        # Debug: Show what's available
        echo "=== Debug Information ==="
        echo "Current working directory: $(pwd)"
        echo "Contents of current directory:"
        ls -la .
        
        echo ""
        echo "=== Checking for mounted directories ==="
        echo "Checking /prsInteractive/data:"
        ls -la /prsInteractive/data/ || echo "❌ /prsInteractive/data not accessible"
        
        echo "Checking /host_prsInteractive (alternative mount):"
        ls -la /host_prsInteractive/data/ || echo "❌ /host_prsInteractive not accessible"
        
        echo ""
        echo "=== Looking for data files in all possible locations ==="
        find / -name "participant.csv" -type f 2>/dev/null | head -5 || echo "No participant.csv found"
        
        # Try multiple data path strategies
        DATA_PATH=""
        
        if [ -f "/prsInteractive/data/participant.csv" ]; then
            echo "✅ Using /prsInteractive/data"
            DATA_PATH="/prsInteractive/data"
        elif [ -f "/host_prsInteractive/data/participant.csv" ]; then
            echo "✅ Using /host_prsInteractive/data"
            DATA_PATH="/host_prsInteractive/data"
        else
            echo "❌ Data files not found in expected locations"
            echo "Available mount points:"
            df -h
            echo ""
            echo "Directory structure from root:"
            find / -maxdepth 3 -type d 2>/dev/null | grep -E "(prs|data)" || echo "No prs or data directories found"
            exit 1
        fi
        
        echo "Using DATA_PATH: $DATA_PATH"
        
        # Verify required files exist
        echo ""
        echo "=== Verifying required files ==="
        for file in participant.csv participant_environment.csv; do
            if [ -f "$DATA_PATH/$file" ]; then
                echo "✅ Found: $file"
                ls -la "$DATA_PATH/$file"
            else
                echo "❌ Missing: $file"
                exit 1
            fi
        done
        
        # Create results directory
        mkdir -p /prsInteractive/results/~{pheno} || mkdir -p ./results/~{pheno}
        RESULTS_PATH="/prsInteractive/results/~{pheno}"
        if [ ! -d "$RESULTS_PATH" ]; then
            RESULTS_PATH="./results/~{pheno}"
        fi
        
        echo "Using RESULTS_PATH: $RESULTS_PATH"
        
        # Find Python script
        echo ""
        echo "=== Finding Python script ==="
        find . -name "create_pheno_train_test_split.py" -type f
        SCRIPT_PATH=$(find . -name "create_pheno_train_test_split.py" -type f | head -1)
        
        if [ -z "$SCRIPT_PATH" ]; then
            echo "❌ Python script not found"
            exit 1
        fi
        
        echo "Using script: $SCRIPT_PATH"
        
        # Run the Python script
        echo ""
        echo "=== Running Python script ==="
        python "$SCRIPT_PATH" \
            --data_path "$DATA_PATH" \
            --icd10 ~{icd10} \
            --pheno ~{pheno} \
            --pheno_str "~{pheno_str}" \
            --pheno_path "$RESULTS_PATH"
            
        # Create completion file
        echo "train test split and phenotype creation complete" > completion_file.txt
        echo "Data path used: $DATA_PATH" >> completion_file.txt
        echo "Results path used: $RESULTS_PATH" >> completion_file.txt
    >>>
    
    runtime {
        docker: "ukb-base:V1"
        memory: "50 GB"
        cpu: 1
    }
    
    output {
        File completion_file = "completion_file.txt"
    }
}