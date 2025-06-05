version 1.0

## Complete prsInteractive Pipeline for DNAnexus
## Replicates the full bash workflow with dynamic data generation

workflow prsInteractivePipeline {
    meta {
        title: "prsInteractive Complete Genomics Pipeline"
        summary: "End-to-end pipeline for UK Biobank phenotype analysis, variant processing, epistasis analysis, and feature selection via ML model "
        description: "This workflow processes UK Biobank data to create phenotype cohorts, clean variant calls, perform epistasis analysis, and train machine learning models for polygenic risk scoring."
        version: "2.0"
        author: "prsInteractive Team"
    }
    
    parameter_meta {
        # Core phenotype parameters
        pheno: "Phenotype name used to create phenotype results (e.g., type2Diabetes, celiacDisease, myocardialInfarction)"
        icd10: "ICD10 code in non-cancer illness code for the phenotype (e.g., E11 for type 2 diabetes, I21 for MI)"
        pheno_str: "String to search for in self-reported illness fields"
        
        # Analysis parameters  
        n_parallel_jobs: "Number of parallel jobs for epistasis analysis (default: 40)"
        platform: "Platform identifier (dnanexus)"
        
        # Input cohort data (dynamic generation instead of pre-uploaded files)
        cohort_participant_data: "UK Biobank participant cohort data"
        cohort_environment_data: "UK Biobank environmental/clinical data"
        cohort_hla_data: "UK Biobank HLA imputation data" 
        cohort_variant_calls: "UK Biobank variant call files (bed/bim/fam format)"
        cohort_withdrawals: "UK Biobank participant withdrawals list"
        
        # Processing scripts
        scripts_archive: "Archive containing all processing scripts"
        
        # Optional parameters
        batch_size: "Batch size for model training (default: 3000)"
        max_model_jobs: "Maximum number of model training jobs to run (default: 2 for testing)"
    }
    
    input {
        # Core parameters
        String pheno = "type2Diabetes"
        String icd10 = "E11"  
        String pheno_str = "type 2 diabetes"
        String platform = "dnanexus"
        Int n_parallel_jobs = 40
        
        # Cohort data files (these replace individual pre-uploaded files)
        File cohort_participant_data
        File cohort_environment_data  
        File cohort_hla_data
        Array[File] cohort_variant_calls
        File cohort_withdrawals
        
        # Scripts
        File scripts_archive
        
        # Optional parameters
        Int batch_size = 3000
        Int max_model_jobs = 2
        
        # Analysis control flags
        Boolean run_epistasis = true
        Boolean run_main_models = true
        Boolean run_epi_models = true
    }
    
    # Step 1: Create phenotype data and train/test/holdout splits
    call createPhenotypeData {
        input:
            pheno = pheno,
            icd10 = icd10,
            pheno_str = pheno_str,
            platform = platform,
            participant_data = cohort_participant_data,
            environment_data = cohort_environment_data,
            scripts_archive = scripts_archive
    }
    
    # Step 2: Clean environmental, HLA, and covariate data
    call cleanEnvironmentalHlaCovarData {
        input:
            pheno = pheno,
            participant_data = cohort_participant_data,
            environment_data = cohort_environment_data,
            hla_data = cohort_hla_data,
            scripts_archive = scripts_archive
    }
    
    # Step 3: Clean variant calls using PLINK
    call cleanVariantCalls {
        input:
            pheno = pheno,
            variant_calls = cohort_variant_calls,
            withdrawals = cohort_withdrawals,
            training_ids = createPhenotypeData.training_ids,
            test_ids = createPhenotypeData.test_ids,
            scripts_archive = scripts_archive
    }
    
    # Step 4: Merge chromosomes and create analysis-ready files
    call mergeChromosomes {
        input:
            pheno = pheno,
            cleaned_variants = cleanVariantCalls.cleaned_variant_files,
            phenotype_labels = createPhenotypeData.phenotype_labels,
            training_ids = createPhenotypeData.training_ids,
            test_ids = createPhenotypeData.test_ids,
            holdout_ids = createPhenotypeData.holdout_ids,
            withdrawals = cohort_withdrawals,
            scripts_archive = scripts_archive
    }
    
    # Step 5: Epistasis analysis (if enabled)
    if (run_epistasis) {
        call multiprocessingFastEpistasis {
            input:
                pheno = pheno,
                n_parallel_jobs = n_parallel_jobs,
                merged_variants = mergeChromosomes.merged_bed_file,
                training_ids = createPhenotypeData.training_ids,
                scripts_archive = scripts_archive
        }
    }
    
    # Step 6: Run main effect models (if enabled)
    if (run_main_models) {
        call runModelBatches as runMainModels {
            input:
                pheno = pheno,
                data_type = "main",
                max_jobs = max_model_jobs,
                batch_size = batch_size,
                training_data = mergeChromosomes.training_raw_file,
                test_data = mergeChromosomes.test_raw_file,
                snp_list = mergeChromosomes.snp_list,
                epi_summary = select_first([multiprocessingFastEpistasis.epi_summary_file, ""]),
                scripts_archive = scripts_archive
        }
    }
    
    # Step 7: Run epistasis models (if enabled and epistasis completed)
    if (run_epi_models && run_epistasis) {
        call runModelBatches as runEpiModels {
            input:
                pheno = pheno,
                data_type = "epi", 
                max_jobs = max_model_jobs,
                batch_size = batch_size,
                training_data = mergeChromosomes.training_raw_file,
                test_data = mergeChromosomes.test_raw_file,
                snp_list = mergeChromosomes.snp_list,
                epi_summary = select_first([multiprocessingFastEpistasis.epi_summary_file]),
                scripts_archive = scripts_archive
        }
    }
    
    output {
        # Core phenotype outputs
        File phenotype_summary = createPhenotypeData.summary_file
        File training_ids = createPhenotypeData.training_ids
        File test_ids = createPhenotypeData.test_ids
        File holdout_ids = createPhenotypeData.holdout_ids
        File phenotype_labels = createPhenotypeData.phenotype_labels
        
        # Processed data outputs
        File cleaned_covar_data = cleanEnvironmentalHlaCovarData.covar_file
        File cleaned_hla_data = cleanEnvironmentalHlaCovarData.hla_file
        File cleaned_environment_data = cleanEnvironmentalHlaCovarData.environment_file
        
        # Variant analysis outputs
        Array[File] cleaned_variant_files = cleanVariantCalls.cleaned_variant_files
        File merged_bed_file = mergeChromosomes.merged_bed_file
        File training_raw_file = mergeChromosomes.training_raw_file
        File test_raw_file = mergeChromosomes.test_raw_file
        File holdout_raw_file = mergeChromosomes.holdout_raw_file
        File snp_list = mergeChromosomes.snp_list
        
        # Epistasis outputs (optional)
        File? epi_summary_file = multiprocessingFastEpistasis.epi_summary_file
        Array[File]? epi_intermediate_files = multiprocessingFastEpistasis.epi_intermediate_files
        
        # Model outputs (optional)
        Array[File]? main_model_scores = runMainModels.model_scores
        Array[File]? main_feature_scores = runMainModels.feature_scores
        Array[File]? main_model_files = runMainModels.model_files
        
        Array[File]? epi_model_scores = runEpiModels.model_scores  
        Array[File]? epi_feature_scores = runEpiModels.feature_scores
        Array[File]? epi_model_files = runEpiModels.model_files
    }
}

task createPhenotypeData {
    meta {
        description: "Creates phenotype data and train/test/holdout splits from UK Biobank participant data"
    }
    
    input {
        String pheno
        String icd10
        String pheno_str
        String platform
        File participant_data
        File environment_data
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Phenotype Data Creation ==="
        echo "Phenotype: ~{pheno}"
        echo "ICD10: ~{icd10}"
        echo "Search string: ~{pheno_str}"
        
        # Create working directories
        mkdir -p data results/~{pheno} logs scripts
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy input data
        cp ~{participant_data} data/participant.csv
        cp ~{environment_data} data/participant_environment.csv
        
        # Activate conda environment
        eval "$(conda shell.bash hook)" || true
        conda activate ukb_env || echo "Warning: Could not activate ukb_env"
        
        # Run phenotype creation script
        python scripts/create_pheno_train_test_split.py \
            --data_path data \
            --pheno_path results/~{pheno} \
            --pheno ~{pheno} \
            --pheno_str "~{pheno_str}" \
            --icd10 ~{icd10}
        
        # Create summary
        cat > workflow_summary.txt << EOF
=== prsInteractive Phenotype Creation Summary ===
Timestamp: $(date)
Phenotype: ~{pheno}
ICD10 Code: ~{icd10}
Search String: ~{pheno_str}
Platform: ~{platform}

Output Files Generated:
$(ls -la results/~{pheno}/)

Processing completed successfully.
EOF
        
        echo "✅ Phenotype data creation completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: 2
        memory: "50 GB"
        disks: "local-disk 20 HDD"
    }
    
    output {
        File summary_file = "workflow_summary.txt"
        File training_ids = "results/~{pheno}/trainingID.txt"
        File test_ids = "results/~{pheno}/testID.txt"
        File holdout_ids = "results/~{pheno}/holdoutID.txt"
        File phenotype_labels = "results/~{pheno}/pheno.txt"
    }
}

task cleanEnvironmentalHlaCovarData {
    meta {
        description: "Cleans environmental, HLA, and covariate data"
    }
    
    input {
        String pheno
        File participant_data
        File environment_data
        File hla_data
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Environmental, HLA, and Covariate Data Cleaning ==="
        
        # Setup directories
        mkdir -p data scripts results
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy input data
        cp ~{participant_data} data/participant.csv
        cp ~{environment_data} data/participant_environment.csv
        cp ~{hla_data} data/hla_participant.csv
        
        # Create HLA headers file (this would normally be provided)
        # For now, create a placeholder - in real implementation, extract from actual data
        echo -e "HLA_A_01:01\tHLA_A_02:01\tHLA_B_07:02" > data/ukb_hla_v2.txt
        
        # Activate conda environment
        eval "$(conda shell.bash hook)" || true
        conda activate ukb_env || echo "Warning: Could not activate ukb_env"
        
        # Run data cleaning script
        python scripts/clean_environment_hla_covar_data.py \
            --data_folder data \
            --results_folder results
        
        echo "✅ Environmental, HLA, and covariate data cleaning completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 30 HDD"
    }
    
    output {
        File covar_file = "results/covar.txt"
        File hla_file = "results/participant_hla.csv"
        File environment_file = "results/participant_environment.csv"
    }
}

task cleanVariantCalls {
    meta {
        description: "Cleans variant calls using PLINK with quality control filters"
    }
    
    input {
        String pheno
        Array[File] variant_calls
        File withdrawals
        File training_ids
        File test_ids
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Variant Call Cleaning ==="
        
        # Setup directories
        mkdir -p data/variant_calls scripts results
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy variant call files
        variant_files=(~{sep=" " variant_calls})
        for file in "${variant_files[@]}"; do
            cp "$file" data/variant_calls/
        done
        
        # Copy other required files
        cp ~{withdrawals} data/withdrawalsID.txt
        cp ~{training_ids} trainingID.txt
        cp ~{test_ids} testID.txt
        
        # Combine training and test IDs
        cat trainingID.txt testID.txt > combinedID.txt
        
        # Set environment variables for script
        export DATA_PATH="$(pwd)/data"
        export PHENO_PATH="$(pwd)/results"
        
        # Create results directory
        mkdir -p results
        
        # Run variant cleaning script
        bash scripts/plink_clean_variant_calls.sh
        
        echo "✅ Variant call cleaning completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: 4
        memory: "32 GB"
        disks: "local-disk 100 HDD"
    }
    
    output {
        Array[File] cleaned_variant_files = glob("results/chr*_cleaned.*")
    }
}

task mergeChromosomes {
    meta {
        description: "Merges chromosome files and creates analysis-ready datasets"
    }
    
    input {
        String pheno
        Array[File] cleaned_variants
        File phenotype_labels
        File training_ids
        File test_ids
        File holdout_ids
        File withdrawals
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Chromosome Merging ==="
        
        # Setup directories
        mkdir -p data scripts results
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy cleaned variant files
        variant_files=(~{sep=" " cleaned_variants})
        for file in "${variant_files[@]}"; do
            cp "$file" results/
        done
        
        # Copy required files
        cp ~{phenotype_labels} results/pheno.txt
        cp ~{training_ids} results/trainingID.txt
        cp ~{test_ids} results/testID.txt
        cp ~{holdout_ids} results/holdoutID.txt
        cp ~{withdrawals} data/withdrawalsID.txt
        
        # Set environment variables
        export DATA_PATH="$(pwd)/data"
        export PHENO_PATH="$(pwd)/results"
        
        # Run chromosome merging script
        bash scripts/merge_chromosomes.sh
        
        echo "✅ Chromosome merging completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: 4
        memory: "64 GB"
        disks: "local-disk 200 HDD"
    }
    
    output {
        File merged_bed_file = "results/merged_allChromosomes.bed"
        File merged_bim_file = "results/merged_allChromosomes.bim"
        File merged_fam_file = "results/merged_allChromosomes.fam"
        File training_raw_file = "results/trainingCombined.raw"
        File test_raw_file = "results/testCombined.raw"
        File holdout_raw_file = "results/holdoutCombined.raw"
        File snp_list = "results/merged_allChromosomes.snplist"
    }
}

task multiprocessingFastEpistasis {
    meta {
        description: "Performs fast epistasis analysis using PLINK with multiprocessing"
    }
    
    input {
        String pheno
        Int n_parallel_jobs
        File merged_bed_file
        File training_ids
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Fast Epistasis Analysis ==="
        echo "Using ~{n_parallel_jobs} parallel jobs"
        
        # Setup directories
        mkdir -p scripts results/epiFiles/preSummaryFiles
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy required files
        cp ~{merged_bed_file} merged_allChromosomes.bed
        cp ~{training_ids} trainingID.txt
        
        # Also need .bim and .fam files (assume they're named consistently)
        cp ${merged_bed_file%.bed}.bim merged_allChromosomes.bim 2>/dev/null || echo "Warning: .bim file not found"
        cp ${merged_bed_file%.bed}.fam merged_allChromosomes.fam 2>/dev/null || echo "Warning: .fam file not found"
        
        # Set environment variables
        export PHENO_PATH="$(pwd)/results"
        export N=~{n_parallel_jobs}
        
        # Run epistasis analysis
        OUT_DIR="results/epiFiles/preSummaryFiles"
        OUTPUT_FILE="${OUT_DIR}/trainingEpi.epi.cc"
        
        echo "Running ~{n_parallel_jobs} parallel PLINK epistasis jobs..."
        
        # Run parallel PLINK jobs
        for i in $(seq 1 ~{n_parallel_jobs}); do
            echo "Starting PLINK job $i"
            plink --bfile merged_allChromosomes \
                  --keep trainingID.txt \
                  --fast-epistasis 'boost' \
                  --parallel $i ~{n_parallel_jobs} \
                  --out "${OUT_DIR}/trainingEpi" &
        done
        
        # Wait for all jobs to complete
        wait
        echo "All parallel PLINK jobs completed"
        
        # Merge results
        > $OUTPUT_FILE
        for i in $(seq 1 ~{n_parallel_jobs}); do
            if [ -f "${OUT_DIR}/trainingEpi.epi.cc.${i}" ]; then
                cat "${OUT_DIR}/trainingEpi.epi.cc.${i}" >> $OUTPUT_FILE
            fi
        done
        
        # Create summary
        plink --epistasis-summary-merge $OUTPUT_FILE ~{n_parallel_jobs} \
              --out "results/epiFiles/trainingCombinedEpi"
        
        echo "✅ Epistasis analysis completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: ~{n_parallel_jobs}
        memory: "200 GB"
        disks: "local-disk 500 HDD"
    }
    
    output {
        File epi_summary_file = "results/epiFiles/trainingCombinedEpi.epi.cc.summary"
        Array[File] epi_intermediate_files = glob("results/epiFiles/preSummaryFiles/*")
    }
}

task runModelBatches {
    meta {
        description: "Runs machine learning model training in batches"
    }
    
    input {
        String pheno
        String data_type  # "main" or "epi"
        Int max_jobs
        Int batch_size
        File training_data
        File test_data
        File snp_list
        String epi_summary
        File scripts_archive
    }
    
    command <<<
        set -euo pipefail
        
        echo "=== Model Batch Training ==="
        echo "Data type: ~{data_type}"
        echo "Max jobs: ~{max_jobs}"
        echo "Batch size: ~{batch_size}"
        
        # Setup directories
        mkdir -p scripts results/scores results/models results/figures
        
        # Extract scripts
        cd scripts
        if [[ "~{scripts_archive}" == *.tar.gz ]] || [[ "~{scripts_archive}" == *.tgz ]]; then
            tar -xzf ~{scripts_archive}
        elif [[ "~{scripts_archive}" == *.zip ]]; then
            unzip ~{scripts_archive}
        fi
        cd ..
        
        # Copy data files
        cp ~{training_data} trainingCombined.raw
        cp ~{test_data} testCombined.raw
        cp ~{snp_list} merged_allChromosomes.snplist
        
        if [ "~{data_type}" == "epi" ] && [ -n "~{epi_summary}" ]; then
            cp ~{epi_summary} epi_summary.txt
            INPUT_FILE="epi_summary.txt"
        else
            INPUT_FILE="merged_allChromosomes.snplist"
        fi
        
        # Calculate number of batches
        TOTAL_LINES=$(wc -l < "$INPUT_FILE")
        TOTAL_BATCHES=$(( (TOTAL_LINES + ~{batch_size} - 1) / ~{batch_size} ))
        
        echo "Total lines: $TOTAL_LINES"
        echo "Total batches: $TOTAL_BATCHES"
        echo "Running first ~{max_jobs} batches for testing"
        
        # Activate conda environment
        eval "$(conda shell.bash hook)" || true
        conda activate ukb_env || echo "Warning: Could not activate ukb_env"
        
        # Set environment variables
        export PHENO_PATH="$(pwd)/results"
        export TRAINING_PATH="$(pwd)/trainingCombined.raw"
        export TEST_PATH="$(pwd)/testCombined.raw"
        export EPI_PATH="$(pwd)/epi_summary.txt"
        export DATA_TYPE="~{data_type}"
        export PHENO="~{pheno}"
        
        # Run model batches (limited for testing)
        BATCHES_TO_RUN=$(( TOTAL_BATCHES < ~{max_jobs} ? TOTAL_BATCHES : ~{max_jobs} ))
        
        for i in $(seq 1 $BATCHES_TO_RUN); do
            echo "Running batch $i"
            
            python scripts/sklearnSectionModelsScoreTrain.py \
                --pheno_folder "results" \
                --training_file "trainingCombined.raw" \
                --test_file "testCombined.raw" \
                --epi_file "epi_summary.txt" \
                --data_type "~{data_type}" \
                --start "$i" \
                --end "$i" \
                --pheno "~{pheno}"
        done
        
        echo "✅ Model batch training completed"
    >>>
    
    runtime {
        docker: "dx://your-project:your-docker-asset"
        cpu: 8
        memory: "64 GB"
        disks: "local-disk 100 HDD"
    }
    
    output {
        Array[File] model_scores = glob("results/scores/*.csv")
        Array[File] feature_scores = glob("results/scores/feature*.csv")
        Array[File] model_files = glob("results/models/*.pkl")
        Array[File] figure_files = glob("results/figures/*")
    }
}