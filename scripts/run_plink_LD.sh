#!/bin/bash

echo "[DIR] scripts directory : $SCRIPTS_DIR"

python "$SCRIPTS_DIR/helper/create_LD_SnpList.py"

wait 

plink --bfile "$PHENO_PATH/merged_allChromosomes" --extract "$PHENO_PATH/finalModelLDSnps.txt" --indep-pairwise 100kb 1 .6 --r2 --show-tags all --out "$PHENO_PATH/finalModel"

#overwites the reducedModelFeatureScores.csv with pruned features in LD and reversed
python "$SCRIPTS_DIR/filter_features_LD.py"