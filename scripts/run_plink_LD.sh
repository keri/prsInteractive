#!/bin/bash

pheno=$1

#PHENO='type2Diabetes_test'


echo "PHENO is set to : $PHENO"



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
	source "${RESULTS_PATH}/$PHENO/pheno.config"
fi

echo "[DIR] scripts directory : $SCRIPTS_DIR"

python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"

wait 

plink --bfile "$PHENO_PATH/merged_allChromosomes" --extract "$PHENO_PATH/finalModelLDSnps.txt" --indep-pairwise 100kb 1 .6 --r2 --show-tags all --out "$PHENO_PATH/finalModel"

#overwites the reducedModelFeatureScores.csv with pruned features in LD and reversed
python "$SCRIPTS_DIR/filter_features_LD.py"