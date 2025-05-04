#!/bin/bash

OUTPUT_FILE="${PHENO_PATH}/mergeChromosomeList.txt"

echo "output file = ${OUTPUT_FILE}..."

# Empty the output file if it already exists
: > $OUTPUT_FILE


# Loop through parallel chunks 1 to 40
for i in $(seq 2 22); 
do
	echo "${PHENO_PATH}/chr${i}_cleaned.bed ${PHENO_PATH}/chr${i}_cleaned.bim ${PHENO_PATH}/chr${i}_cleaned.fam" >> $OUTPUT_FILE

done

plink --bfile "${PHENO_PATH}"/chr1_cleaned --merge-list "${OUTPUT_FILE}" --write-snplist --make-bed --out "${PHENO_PATH}"/merged_allChromosomes

# Empty the output file if it already exists
: > $OUTPUT_FILE

# Loop through parallel chunks 1 to 40
for i in $(seq 2 22); 
	
do
	echo "${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.bed ${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.bim ${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.fam" >> $OUTPUT_FILE
	
done

plink --bfile "${DATA_PATH}"/variant_calls/ukb22418_c1_b0_v2 --merge-list "${OUTPUT_FILE}" --keep "${PHENO_PATH}"/holdoutID.txt --pheno "${PHENO_PATH}"/pheno.txt --extract "${PHENO_PATH}"/merged_allChromosomes.snplist --make-bed --out "${PHENO_PATH}"/holdoutCombined

rm $OUTPUT_FILE

# Loop through parallel chunks 1 to 40
for i in $(seq 1 22); 
	
do

	rm "${PHENO_PATH}/chr${i}_cleaned"*
	
done
