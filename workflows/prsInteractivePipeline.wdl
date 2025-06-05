version 1.0

workflow prsInteractivePipeline {
    input {
        String pheno
        String icd10
        String pheno_str
        Int n
        Array[File] raw_variant_calls
        File participant_hla
        File participant_data
        File participant_environment
        File covar
        File hla_headers
        File withdrawals
        File config_file
    }
    
    # Set up the environment variables and create pheno folder
    call createSummary {
        input:
            pheno = pheno,
            icd10 = icd10,
            pheno_str = pheno_str,
            config_file = config_file
    }
    
}

task createSummary { 
    input {
        String pheno
        String icd10
        String pheno_str
        File config_file
        
    }
    
    command <<<
        set -e
        
        echo "Processing phenotype: ~{pheno}"
        echo "Using pheno string: ~{pheno_str}"
        echo "With icd10: ~{icd10}"
        
        # Debug information
        echo "Current working directory: $(pwd)"
        echo "Current user: $(whoami)"
        ls -la . || true
        
        # Create the summary file to capture time stamp, phenotype substring, and icd10 code
        ts=$(date +%Y%m%d_%H%M%S)
        
        # Create summary file in working directory
        summaryFile="summary.txt"
        touch "${summaryFile}"
        echo "timestamp of analysis: ${ts}" >> "${summaryFile}"
        echo "phenotype analyzed: ~{pheno}" >> "${summaryFile}"
        echo "phenotype substring filtered in participant.csv: ~{pheno_str}" >> "${summaryFile}"
        echo "icd 10 code used to define phenotype: ~{icd10}" >> "${summaryFile}"
        
        # Verify file was created
        echo "Summary file contents:"
        cat "${summaryFile}"
    >>>

    runtime {
        docker: "ukb-base:V1"
    }
    
    output {
        File summary_file = "summary.txt"
    }

}