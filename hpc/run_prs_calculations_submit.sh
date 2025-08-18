#!/bin/bash -x

#
#SBATCH --job-name=calculate_prs
#SBATCH -o /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_prs.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_calculate_prs.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=oconnoas@staff.vuw.ac.nz
#



module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env

pheno=$1
#pheno="celiacDisease"


# Source config with error handling
if [ ! -f "../env.config" ]; then
    echo "ERROR: ../env.config not found!"
    echo "Current directory: $(pwd)"
    echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
    exit 1
    
else
    source ../env.config
fi

#check that a results folder for phenotype exists
if [ ! -d "${RESULTS_PATH}/$pheno" ]; then
    echo "Folder '${RESULTS_PATH}/$pheno' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${RESULTS_PATH}/$pheno/pheno.config"
fi



# Check if RESULTS_PATH was set from env.config
if [ -z "$RESULTS_PATH" ]; then
    echo "ERROR: RESULTS_PATH not set from env.config"
    exit 1
fi

echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"


# Set phenotype-specific paths


echo "[WORKFLOW] ENV_TYPE is set to: $env_type"
export PHENO_PATH="$RESULTS_PATH/$pheno"
export TEST_PATH=$TEST_PATH
export HOLDOUT_PATH=$HOLDOUT_PATH
export PHENO=$pheno
export RESULTS_PATH=$RESULTS_PATH
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH
export DATA_PATH=$DATA_PATH
export HLA_FILE=$HLA_FILE
export COVAR_FILE=$COVAR_FILE
export GENE_ENV_TEST=$GENE_ENV_TEST
export GENE_ENV_HOLDOUT=$GENE_ENV_HOLDOUT
export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
export ENV_FILE=$ENV_FILE


echo "[DEBUG] ===== ENVIRONMENT VARIABLES ====="
echo "PHENOTYPE BEING ANALYZED: $PHENO"
echo "PHENO_PATH: $PHENO_PATH"
echo "SCRIPTS_DIR: $SCRIPTS_DIR"
echo "TEST_PATH: $TEST_PATH"
echo "HOLDOUT_PATH: $HOLDOUT_PATH"
echo "RESULTS_PATH: $RESULTS_PATH"
echo "TEST_ENV_GEN_FILE: $GENE_ENV_TEST"
echo "HOLDOUT_ENV_GEN_FILE: $GENE_ENV_HOLDOUT"
echo "====================================="

# Check if required files exist before running Python
echo "[DEBUG] Checking required files..."

required_files=(
    "$FEATURE_SCORES_FILE"
    "$GENE_ENV_HOLDOUT"
    "$GENE_ENV_TEST"
    "$WITHDRAWAL_PATH"
    "$COVAR_FILE"
    "$HLA_FILE"
    "$HOLDOUT_PATH"
    "$TEST_PATH"
    "${SCRIPTS_DIR}/calculate_prs_for_filtered_main_epi.py"
    "$SCRIPTS_DIR/run_plink_LD.sh"
    "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file does not exist: $file"
        exit 1
    else
        echo "[DEBUG] Found: $file"
    fi
done

echo "[DEBUG] All required files found. Starting Python script..."

export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
#check to see if LD has been done previously before association
if [ ! -f "$PHENO_PATH/finalModel.ld" ];then
    if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
        #run LD script
        export PRE_POST_ASSOCIATION='post'
    else
        export PRE_POST_ASSOCIATION='pre'
    fi
    sbatch run_plink_LD_submit.sh $pheno
    sleep 15m # do nothing for 15 minutes while LD is done 
fi

#check to see if gene-environment additive analysis has been done
if [ ! -f "$PHENO_PATH/scores/featureScoresReducedFinalModel.filtered.csv" ]; then
    
    export SCORES_PATH="${PHENO_PATH}/scores"
    python "$SCRIPTS_DIR/filter_non_additive_gen_env_features.py"
    
    #ensure file is there
    if [ ! -f "$PHENO_PATH/scores/featureScoresReducedFinalModel.filtered.csv" ]; then
        echo "❌ Expected filtered file not found: $PHENO_PATH/scores/featureScoresReducedFinalModel.filtered.csv"
        echo "Continuing with original FEATURE_SCORES_FILE"
    else
        #update config file
        NEW_FEATURE_SCORES_FILE="$PHENO_PATH/scores/featureScoresReducedFinalModel.filtered.csv" 
        #check to see if config has been updated
        if [ "$FEATURE_SCORES_FILE" == "featureScoresReducedFinalModel.csv" ]; then
            # Update config file
            CONFIG_FILE="${PHENO_PATH}/pheno.config"
            echo "Updating config file: $CONFIG_FILE"
            
            # Create backup
            cp "$CONFIG_FILE" "${CONFIG_FILE}.backup"
            
            if [[ "$(uname)" == "Darwin" ]]; then
                # macOS version
                sed -i '' "s|^FEATURE_SCORES_FILE=.*|FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}|" "$CONFIG_FILE"
            else
                # Linux version  
                sed -i "s|^FEATURE_SCORES_FILE=.*|FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}|" "$CONFIG_FILE"
            fi
            
            # Verify the change
            if grep -q "^FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}$" "$CONFIG_FILE"; then
                echo "✓ Successfully updated EPI_FILE in config"
                rm "${CONFIG_FILE}.backup"  # Remove backup if successful
                FEATURE_SCORES_FILE=$EW_FEATURE_SCORES_FILE
                export FEATURE_SCORES_FILE=$NEW_FEATURE_SCORES_FILE
            else
                echo "❌ Failed to update config file"
                mv "${CONFIG_FILE}.backup" "$CONFIG_FILE"  # Restore backup
                
            fi
        fi
    fi
fi




# Run the Python script
python "${SCRIPTS_DIR}/calculate_prs_for_filtered_main_epi.py"

exit_code=$?
echo "[DEBUG] Python script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python script failed with exit code $exit_code"
    exit $exit_code
fi

python "${SCRIPTS_DIR}/combine_prs.py"

exit_code=$?
echo "[DEBUG] Python combine_prs.py script exited with code: $exit_code"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python combine_prs.py script failed with exit code $exit_code"
    exit $exit_code
fi

#need exported PHENO_PATH
#calculate performance metrics for models using yProba and prs calculations for validation and holdout data
python "${SCRIPTS_DIR}/calculate_prs_stats.py"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Python calculate_prs_stats.py script failed with exit code $exit_code"
    exit $exit_code
fi

echo "[DEBUG] Script completed successfully"












