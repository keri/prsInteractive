#!/bin/bash

pheno=$1

#pheno='type2Diabetes'


# Source config with error handling
if [ ! -f "../env.config" ]; then
	echo "ERROR: ../env.config not found!"
	echo "Current directory: $(pwd)"
	echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
	exit 1
	
else
	source ../env.config
fi

PHENO_PATH="${RESULTS_PATH}/$pheno"

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

echo "[DIR] scripts directory : $SCRIPTS_DIR"


export PHENO_PATH=$PHENO_PATH
python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"

wait 

plink --bfile "$PHENO_PATH/merged_allChromosomes" --extract "$PHENO_PATH/finalModelLDSnps.txt" --indep-pairwise 100kb 1 .6 --r2 --show-tags all --out "$PHENO_PATH/finalModel"

#creates new feature_scores_file with pruned features in LD and reversed
export PRE_POST_ASSOCIATION=$PRE_POST_ASSOCIATION
export FEATURE_SCORES_FILE=$FEATURE_SCORES_FILE
python "$SCRIPTS_DIR/filter_features_LD.py"

##update FEATURE_SCORES_FILE conditional on previous steps: featureScoresReducedFinalModel.csv post association analysis else before
#if [ "$FEATURE_SCORES_FILE" == "$PHENO_PATH/scores/featureScoresReducedFinalModel.csv" ]; then
#	NEW_FEATURE_SCORES_FILE="$FEATURE_SCORES_FILE"
#else
#	NEW_FEATURE_SCORES_FILE="$PHENO_PATH/scores/importantFeaturesForAssociationAnalysis.csv"
#fi
## Verify filtered file exists
#if [ ! -f "$NEW_FEATURE_SCORES_FILE" ]; then
#	echo "❌ Expected filtered file not found: $NEW_FEATURE_SCORES_FILE"
#	echo "Continuing with original FEATURE_SCORES_FILE file"
#
#else
#	# Update config file
#	CONFIG_FILE="${PHENO_PATH}/pheno.config"
#	echo "Updating config file: $CONFIG_FILE"
#	echo "New FEATURE_SCORES_FILE value: $NEW_FEATURE_SCORES_FILE"
#	
#	# Create backup
#	cp "$CONFIG_FILE" "${CONFIG_FILE}.backup"
#	
#	if [[ "$(uname)" == "Darwin" ]]; then
#		# macOS version
#		sed -i '' "s|^FEATURE_SCORES_FILE=.*|FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}|" "$CONFIG_FILE"
#	else
#		# Linux version  
#		sed -i "s|^FEATURE_SCORES_FILE=.*|FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}|" "$CONFIG_FILE"
#	fi
#	
#	# Verify the change
#	if grep -q "^FEATURE_SCORES_FILE=${NEW_FEATURE_SCORES_FILE}$" "$CONFIG_FILE"; then
#		echo "✓ Successfully updated FEATURE_SCORES_FILE in config"
#		rm "${CONFIG_FILE}.backup"  # Remove backup if successful
#	else
#		echo "❌ Failed to update config file"
#		mv "${CONFIG_FILE}.backup" "$CONFIG_FILE"  # Restore backup
#		
#	fi
#fi
