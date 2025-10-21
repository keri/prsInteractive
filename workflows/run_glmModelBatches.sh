#!/usr/bin/env

pheno=$1
#pheno='myocardialInfarction'
threshold=${2:-2}  # Default threshold of 2 if not provided
#threshold=4

echo "Using z-score threshold: $threshold"

# Source config
source ../env.config # because you're in prsInteractive/hpc


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

export RESULTS_PATH=$RESULTS_PATH
export PHENO_PATH=$PHENO_PATH
export DATA_PATH=$DATA_PATH
export HLA_FILE="$RESULTS_PATH/participant_hla.csv"
export COVAR_FILE="$RESULTS_PATH/covar.csv"
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export TRAINING_ENV_GEN_FILE=$GENE_ENV_TRAINING
export TEST_ENV_GEN_FILE=$GENE_ENV_TEST

### CHECK TO SEE IF LD HAS BEEN DONE

if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
    echo "LD not completed for importantFeaturesPostShap.csv"
    #run LD script
    export PRE_POST_ASSOCIATION='pre'
    bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
fi


#check to see LD was completed
if [ ! -f "$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv" ]; then
    echo "LD analysis failed for importantFeaturesPostShap.csv"
    FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesPostShap.csv"
    echo "reverting to original importantFeaturesPostShap.csv"
else
    echo "LD analysis completed for importantFeaturesPostShap.csv"
    FEATURES_FOR_MODEL_FILE="$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv"
fi

### NEW PREFILTERING STEP - Filter features with z_score > 3
echo "Starting feature prefiltering step..."

# Check if the features file has a z_score column
if [ -f "$FEATURES_FOR_MODEL_FILE" ]; then
    # Check if shap_zscore column exists
    if head -1 "$FEATURES_FOR_MODEL_FILE" | grep -q "shap_zscore"; then
        echo "shap_zscore column found. Filtering features with |shap_zscore| > $threshold..."
        
        # Create filtered file name
        FILTERED_FEATURES_FILE="$PHENO_PATH/scores/importantFeaturesFiltered_zscore${threshold}.csv"
        
        # First, find the shap_zscore column index
        Z_COL=$(head -1 "$FEATURES_FOR_MODEL_FILE" | tr ',' '\n' | grep -n "shap_zscore" | cut -d: -f1 | head -1)
        
        if [ -z "$Z_COL" ]; then
            echo "Error: shap_zscore column not found"
            echo "Available columns:"
            head -1 "$FEATURES_FOR_MODEL_FILE" | tr ',' '\n' | nl
            exit 1
        fi
        
        echo "Found shap_zscore column at position: $Z_COL"
        
        # Filter using awk with the found column index
        awk -F',' -v zcol="$Z_COL" -v thresh="$threshold" '
        NR==1 {print; next} 
        NR>1 {
            z_val = $zcol;
            gsub(/["'"'"']/, "", z_val);
            if (z_val != "" && (z_val > thresh || z_val < -thresh)) {
                print;
            }
        }' "$FEATURES_FOR_MODEL_FILE" > "$FILTERED_FEATURES_FILE"
        
        # Check if filtering was successful
        if [ $? -eq 0 ] && [ -f "$FILTERED_FEATURES_FILE" ]; then
            # Count original vs filtered features
            ORIGINAL_COUNT=$(tail -n +2 "$FEATURES_FOR_MODEL_FILE" | wc -l)
            FILTERED_COUNT=$(tail -n +2 "$FILTERED_FEATURES_FILE" | wc -l)
            
            echo "Original features: $ORIGINAL_COUNT"
            echo "Filtered features (|z_score| > $threshold): $FILTERED_COUNT"
            
            if [ "$FILTERED_COUNT" -gt 0 ]; then
                echo "Using filtered features file for modeling"
                FEATURES_FOR_MODEL_FILE="$FILTERED_FEATURES_FILE"
            else
                echo "Warning: No features pass z_score > $threshold filter. Using original file."
            fi
        else
            echo "Warning: Filtering failed. Using original features file."
        fi
    else
        echo "No z_score column found in $FEATURES_FOR_MODEL_FILE"
        echo "Available columns:"
        head -1 "$FEATURES_FOR_MODEL_FILE"
        echo "Proceeding with original file..."
    fi
else
    echo "Features file not found: $FEATURES_FOR_MODEL_FILE"
    exit 1
fi

echo "Feature filtering completed. Using: $FEATURES_FOR_MODEL_FILE"

#bash "$SCRIPTS_DIR/run_plink_LD.sh" $pheno
echo "Starting R script..."

# Pass to R script
Rscript "$SCRIPTS_DIR/glmFinalModellingBatches.R" \
--results_path "$RESULTS_PATH" \
--data_path "$DATA_PATH" \
--hla_file "$HLA_FILE" \
--covar_file "$COVAR_FILE" \
--pheno_path "$PHENO_PATH" \
--training_file "$TRAINING_PATH" \
--test_file "$TEST_PATH" \
--training_env_gen_file "$GENE_ENV_TRAINING" \
--test_env_gen_file "$GENE_ENV_TEST" \
--feature_model_file "$FEATURES_FOR_MODEL_FILE"

R_EXIT_CODE=$?

if [ $R_EXIT_CODE -ne 0 ]; then
    echo "ERROR: R script failed with exit code $R_EXIT_CODE"
    exit 1
fi

echo "R script completed successfully. Checking outputs..."

# Check if R script outputs exist before updating config
missing_files=()

if [ ! -f "$PHENO_PATH/scores/featureScoresReducedFinalModelBatches.csv" ]; then
    missing_files+=("featureScoresReducedFinalModelBatches.csv")
fi

if [ ! -f "$PHENO_PATH/scores/modelScoresReducedFinalModelBatches.csv" ]; then
    missing_files+=("modelScoresReducedFinalModelBatches.csv")
fi

if [ ! -f "$PHENO_PATH/scores/predictProbsReducedFinalModelBatches.csv" ]; then
    missing_files+=("predictProbsReducedFinalModelBatches.csv")
fi

# Only update config if all files exist
if [ ${#missing_files[@]} -eq 0 ]; then
    echo "All R script outputs found. Updating pheno.config..."
    {
        echo "FEATURE_SCORES_FILE=$PHENO_PATH/scores/featureScoresReducedFinalModelBatches.csv"
        echo "FINAL_MODEL_SCORES=$PHENO_PATH/scores/modelScoresReducedFinalModelBatches.csv"
        echo "FINAL_MODEL_PROBABILITIES=$PHENO_PATH/scores/predictProbsReducedFinalModelBatches.csv"
    } >> "${RESULTS_PATH}/$pheno/pheno.config"
else
    echo "ERROR: R script failed to generate the following files:"
    printf '%s\n' "${missing_files[@]}"
    echo "pheno.config not updated"
    exit 1
fi

sbatch "run_prs_calculations_submit.sh" $pheno