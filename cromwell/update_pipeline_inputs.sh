#!/bin/bash

# Dynamic Pipeline Input Updater
# Usage: ./update_pipeline_inputs.sh <phenotype> <icd10> <pheno_string> <n_value>

set -e

# Function to display usage
usage() {
    echo "Usage: $0 <phenotype> <icd10> <pheno_string> <n_value> <platform default=local>"
    echo "Example: $0 type2Diabetes E11 'type 2 diabetes' 40 hpc"
    exit 1
}

# Check if correct number of arguments provided
if [ $# -ne 5 ]; then
    echo "Error: Incorrect number of arguments"
    usage
fi

# Assign arguments to variables
PHENOTYPE=$1
ICD10=$2
PHENO_STRING=$3
N_VALUE=$4
PLATFORM=$5


# Get the project root directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_INPUTS_FILE="$PROJECT_ROOT/workflows/pipelineInputs.json"


#create pheno config and .env files for workflow
PHENO_DIR="$PROJECT_ROOT/results/$PHENOTYPE"

echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENOTYPE"
echo "PHENOTYPE STRING USED IN ANALYSIS ...: $PHENO_STRING"
echo "ICD10 USED IN ANALYSIS ...: $ICD10"
echo "PLATFORM USED ON BACKEND ..: $PLATFORM"

#check to see if PHENO_PATH is already present
#create a phenotype env for later use
# Create a .env file for future reference
cat > "${PHENO_DIR}/pheno.config" << EOF
#pheno environment variables
PHENO="${pheno}"
PHENO_PATH=$PHENO_DIR
PHENO_STR="${phenoStr}"
ICD10=$icd10
N=$n
EOF

echo "Updating pipeline inputs for phenotype: $PHENOTYPE"
echo "ICD10 code: $ICD10"
echo "Phenotype string: $PHENO_STRING"
echo "N value: $N_VALUE"

# Create backup of original file
if [ -f "$PIPELINE_INPUTS_FILE" ]; then
    cp "$PIPELINE_INPUTS_FILE" "${PIPELINE_INPUTS_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
    echo "Backup created: ${PIPELINE_INPUTS_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
fi

# Generate the updated pipelineInputs.json
cat > "$PIPELINE_INPUTS_FILE" << EOF
{
    "prsInteractivePipeline.pheno": "$PHENOTYPE",
    "prsInteractivePipeline.icd10": "$ICD10",
    "prsInteractivePipeline.pheno_str": "$PHENO_STRING",
    "prsInteractivePipeline.n": $N_VALUE,
    "prsInteractivePipeline.scripts_archive": "$PROJECT_ROOT/scripts.tar.gz",
    "prsInteractivePipeline.root_directory": "$PROJECT_ROOT",
    "prsInteractivePipeline.platform": "$PLATFORM",
    "prsInteractivePipeline.raw_variant_calls": [
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c1_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c2_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c3_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c4_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c5_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c6_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c7_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c8_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c9_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c10_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c11_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c12_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c13_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c14_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c15_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c16_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c17_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c18_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c19_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c20_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c21_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c22_b0_v2",
        "$PROJECT_ROOT/data/variant_calls/ukb22418_c23_b0_v2"
    ],
    "prsInteractivePipeline.participant_data": "$PROJECT_ROOT/data/participant.csv",
    "prsInteractivePipeline.participant_environment": "$PROJECT_ROOT/data/participant_environment.csv",
    "prsInteractivePipeline.covar": "$PROJECT_ROOT/data/covar.txt",
    "prsInteractivePipeline.hla_headers": "$PROJECT_ROOT/data/ukb_hla_v2.txt",
    "prsInteractivePipeline.participant_hla": "$PROJECT_ROOT/data/hla_participant.csv",
    "prsInteractivePipeline.withdrawals": "$PROJECT_ROOT/data/withdrawals.csv",
    "prsInteractivePipeline.config_file": "$PROJECT_ROOT/config/$PLATFORM.config"
}
EOF

echo "✓ pipelineInputs.json updated successfully"

## Now run the environment setup
#echo "Running environment setup..."
#if [ -f "$PROJECT_ROOT/envSetUp.sh" ]; then
#   chmod +x "$PROJECT_ROOT/envSetUp.sh"
#   "$PROJECT_ROOT/envSetUp.sh" "$PHENOTYPE" "$ICD10" "$PHENO_STRING" "$N_VALUE"
#   echo "✓ Environment setup completed"
#else
#   echo "Warning: envSetUp.sh not found in $PROJECT_ROOT"
#fi
rm "${PIPELINE_INPUTS_FILE}.backup"*
echo "Pipeline configuration updated for phenotype: $PHENOTYPE"