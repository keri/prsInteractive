#!/bin/bash

pheno=$1
epi_subfolder=${2:-""}   # optional: "productEpi" or "summedEpi"
                          # omit to run LD at the pheno root (default: post-epi discovery, pre gene-environment)

# Source global config
if [ ! -f "../env.config" ]; then
	echo "ERROR: ../env.config not found!"
	echo "Current directory: $(pwd)"
	echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
	exit 1
fi
source ../env.config

# Construct and validate PHENO_PATH
PHENO_PATH="${RESULTS_PATH}/${pheno}"
if [ ! -d "${PHENO_PATH}" ]; then
	echo "Folder '${PHENO_PATH}' does not exist..."
	echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
	exit 1
fi

# Default: PHENO_DATA = pheno root (LD run before gene-environment feature discovery)
# With epi_subfolder: PHENO_DATA = subfolder, sources its pheno.config for additional env vars
if [ -n "${epi_subfolder}" ]; then
	PHENO_DATA="${PHENO_PATH}/${epi_subfolder}"
	if [ ! -d "${PHENO_DATA}" ]; then
		echo "Subfolder '${PHENO_DATA}' does not exist..."
		exit 1
	fi
	source "${PHENO_DATA}/pheno.config"
else
	PHENO_DATA="${PHENO_PATH}"
fi

echo "PHENO_PATH: ${PHENO_PATH}"
echo "PHENO_DATA: ${PHENO_DATA}"
echo "[DIR] scripts directory : $SCRIPTS_DIR"

# Build SNP list from importantFeaturesPostShap.csv
export PHENO_DATA
python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"

wait

# Compute pairwise LD — genotype data always at pheno root, output in PHENO_DATA
plink --bfile "$PHENO_PATH/merged_allChromosomes" \
      --extract "$PHENO_DATA/finalModelLDSnps.txt" \
      --indep-pairwise 100kb 1 .6 --r2 --show-tags all \
      --out "$PHENO_DATA/finalModel"

# Filter features for LD
FEATURES_TO_FILTER_LD="${PHENO_DATA}/scores/importantFeaturesPostShap.csv"
if [ ! -f "${FEATURES_TO_FILTER_LD}" ]; then
	echo "Feature file for LD filtering not found: ${FEATURES_TO_FILTER_LD}"
	exit 1
fi

export PHENO_PATH
export PHENO_DATA
export FEATURES_TO_FILTER_LD
python "$SCRIPTS_DIR/filter_features_LD.py"
