#!/bin/bash

cat "${PHENO_PATH}"/testID.txt "${PHENO_PATH}"/trainingID.txt > "${PHENO_PATH}"/combinedID.txt

#chr6 is treated differently as the MHC region is filtered out
plink --bfile  "${DATA_PATH}"/variant_calls/ukb22418_c6_b0_v2 --pheno "${PHENO_PATH}"/pheno.txt --exclude-snps null_26-null_28 --keep "${PHENO_PATH}"/combinedID.txt --remove "${DATA_PATH}"/withdrawalsID.txt --geno .05 --maf .0001 --hwe .00005  --make-bed --out  "${PHENO_PATH}"/chr6_cleaned

for ((a=1; a<23; a++)); 
do
	plink --bfile  "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2 --pheno "${PHENO_PATH}"/pheno.txt --keep "${PHENO_PATH}"/combinedID.txt --remove "${DATA_PATH}"/withdrawalsID.txt --geno .05 --maf .0001 --hwe .00005  --make-bed --out  "${PHENO_PATH}"/chr"${a}"_cleaned

done

for ((a=7; a<23; a++));
do
	plink --bfile  "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2 --pheno "${PHENO_PATH}"/pheno.txt --keep "${PHENO_PATH}"/combinedID.txt --remove "${DATA_PATH}"/withdrawalsID.txt --geno .05 --maf .0001 --hwe .00005 --make-bed --out  "${PHENO_PATH}"/chr"${a}"_cleaned
done



