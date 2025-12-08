#!/bin/bash

pheno=$1

#pheno='type2Diabetes'

#PRE_POST_ASSOCIATION=${2:-"pre"}



# Source config with error handling
if [ ! -f "../env.config" ]; then
	echo "ERROR: ../env.config not found!"
	echo "Current directory: $(pwd)"
	echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
	exit 1
	
else
	source ../env.config
fi

echo "PHENO_DATA: ${PHENO_DATA}"


#FEATURE_SCORES_FILE=${FEATURE_SCORES_FILE:-"${PHENO_PATH}/scores/importantFeaturesPostShap.csv"}

#check that a results folder for phenotype exists
if [ ! -d "${PHENO_PATH}" ]; then
	echo "Folder '${PHENO_PATH}' does not exist..."
	echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
	exit 1
	
else
	echo "sourcing $pheno env variables."
	#source pheno specific environment variables
	source "${PHENO_DATA}/pheno.config"
fi

echo "[DIR] scripts directory : $SCRIPTS_DIR"


export PHENO_DATA
python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"

wait 

#creates the LD SNP data files to be used after looking for important features step
plink --bfile "$PHENO_PATH/merged_allChromosomes" --extract "$PHENO_DATA/finalModelLDSnps.txt" --indep-pairwise 100kb 1 .6 --r2 --show-tags all --out "$PHENO_DATA/finalModel"

#creates new feature_scores_file with pruned features in LD and reversed

export PHENO_PATH
FEATURES_TO_FILTER_LD="${PHENO_DATA}/scores/importantFeaturesPostShap.csv"
#ensure the reduced features file is present
if [ ! -f "${FEATURES_TO_FILTER_LD}" ]; then
	echo "reduced feature file to use in LD analysis is not present"
	exit 1

fi


export FEATURES_TO_FILTER_LD
python "$SCRIPTS_DIR/filter_features_LD.py"


