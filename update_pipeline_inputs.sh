#!/bin/bash

# Dynamic Pipeline Input Updater
# Usage: ./update_pipeline_inputs.sh <phenotype> <icd10> <pheno_string> <n_value>

set -e

# Function to display usage
usage() {
    echo "Usage: $0 <phenotype> <icd10> <pheno_string> <n_value>"
    echo "Example: $0 type2Diabetes E11 'type 2 diabetes' 40"
    exit 1
}

# Check if correct number of arguments provided
if [ $# -ne 4 ]; then
    echo "Error: Incorrect number of arguments"
    usage
fi

# Assign arguments to variables
PHENOTYPE=$1
ICD10=$2
PHENO_STRING=$3
N_VALUE=$4

# Get the project root directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_INPUTS_FILE="$PROJECT_ROOT/pipelineInputs.json"

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
    "prsInteractivePipeline.raw_variant_calls": [
        "/prsInteractive/data/variant_calls/ukb22418_c1_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c2_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c3_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c4_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c5_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c6_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c7_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c8_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c9_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c10_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c11_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c12_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c13_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c14_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c15_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c16_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c17_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c18_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c19_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c20_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c21_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c22_b0_v2",
        "/prsInteractive/data/variant_calls/ukb22418_c23_b0_v2"
    ],
    "prsInteractivePipeline.hla_participant_csv": "/prsInteractive/data/hla_participant.csv",
    "prsInteractivePipeline.participant_data": "/prsInteractive/data/participant.csv",
    "prsInteractivePipeline.participant_environment": "/prsInteractive/data/participant_environment.csv",
    "prsInteractivePipeline.covar": "/prsInteractive/data/covar.txt",
    "prsInteractivePipeline.hla_headers": "/prsInteractive/data/ukb_hla_v2.txt",
    "prsInteractivePipeline.participant_hla": "/prsInteractive/data/hla_participant.csv",
    "prsInteractivePipeline.withdrawals": "/prsInteractive/data/withdrawals.csv",
    "prsInteractivePipeline.config_file": "config/default.config"
}
EOF

echo "✓ pipelineInputs.json updated successfully"

# Now run the environment setup
echo "Running environment setup..."
if [ -f "$PROJECT_ROOT/envSetUp.sh" ]; then
    chmod +x "$PROJECT_ROOT/envSetUp.sh"
    "$PROJECT_ROOT/envSetUp.sh" "$PHENOTYPE" "$ICD10" "$PHENO_STRING" "$N_VALUE"
    echo "✓ Environment setup completed"
else
    echo "Warning: envSetUp.sh not found in $PROJECT_ROOT"
fi

echo "Pipeline configuration updated for phenotype: $PHENOTYPE"