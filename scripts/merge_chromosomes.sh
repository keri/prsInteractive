#!/bin/bash


plink --bfile "${DATA_PATH}"/chr1_cleaned --merge-list "${DATA_PATH}"/mergeChromosomes.txt  --write-snplist --make-bed --out "${DATA_PATH}"/merged_allChromosomes

plink --bfile "${DATA_PATH}"/variant_calls/ukb22418_c1_b0_v2 --merge-list "${DATA_PATH}"/variant_calls/mergeChomosomesHoldout.txt --keep "${PHENO_PATH}"/holdoutID.txt --pheno "${PHENO_PATH}"/pheno.txt --extract "${DATA_PATH}"/merged_allChromosomes.snplist --make-bed --out "${PHENO_PATH}"/holdoutCombined
