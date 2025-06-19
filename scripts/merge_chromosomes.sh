#!/bin/bash

OUTPUT_FILE="${PHENO_PATH}/mergeChromosomeList.txt"

echo "output file = ${OUTPUT_FILE}..."

# Empty the output file if it already exists
: > $OUTPUT_FILE


# Loop through parallel chunks 1 to 22 for combined training test set
for i in $(seq 2 22); 
do
	echo "${PHENO_PATH}/chr${i}_cleaned.bed ${PHENO_PATH}/chr${i}_cleaned.bim ${PHENO_PATH}/chr${i}_cleaned.fam" >> $OUTPUT_FILE

done

plink --bfile "${PHENO_PATH}"/chr1_cleaned --merge-list "${OUTPUT_FILE}" --pheno "${PHENO_PATH}"/pheno.txt --write-snplist --make-bed --out "${PHENO_PATH}"/merged_allChromosomes

# Empty the output file if it already exists
: > $OUTPUT_FILE

# Loop through parallel chunks 1 to 22 for holdout se
for i in $(seq 2 22); 
	
do
	echo "${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.bed ${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.bim ${DATA_PATH}/variant_calls/ukb22418_c${i}_b0_v2.fam" >> $OUTPUT_FILE
	
done

plink --bfile "${DATA_PATH}"/variant_calls/ukb22418_c1_b0_v2 --merge-list "${OUTPUT_FILE}" --keep "${PHENO_PATH}"/holdoutID.txt --remove "${DATA_PATH}"/withdrawalsID.txt --pheno "${PHENO_PATH}"/pheno.txt --extract "${PHENO_PATH}"/merged_allChromosomes.snplist --make-bed --out "${PHENO_PATH}"/holdoutCombined

rm $OUTPUT_FILE

# Loop through parallel chunks 1 to 40
for i in $(seq 1 22); 
	
do

	rm "${PHENO_PATH}/chr${i}_cleaned"*
	
done




####################### CONVERT TO .RAW FILES TO BE USED IN MODELLING  #############################

plink --bfile "${PHENO_PATH}"/holdoutCombined --recode A --keep "${PHENO_PATH}"/holdoutID.txt --out "${PHENO_PATH}"/holdoutCombinedRaw


plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}"/trainingID.txt --out "${PHENO_PATH}"/trainingCombinedRaw

plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}"/testID.txt --out "${PHENO_PATH}"/testCombinedRaw


####################### CONVERT COLUMNS IN FILES TO REMOVE THE VARIANTS IN ORDER TO MATCH OUTPUT FROM EPISTATIC ANALYSIS  #############################

awk 'NR==1 {
	for (i=1; i<=NF; i++) {
		sub(/_[^_]*$/, "", $i)
	}
} { print }' OFS='\t' "${PHENO_PATH}/holdoutCombinedRaw.raw" > "${PHENO_PATH}/holdoutCombined.raw"
	
rm "${PHENO_PATH}/holdoutCombinedRaw.raw"



awk 'NR==1 {
	for (i=1; i<=NF; i++) {
		sub(/_[^_]*$/, "", $i)
	}
} { print }' OFS='\t' "${PHENO_PATH}/trainingCombinedRaw.raw" > "${PHENO_PATH}/trainingCombined.raw"


rm "${PHENO_PATH}/trainingCombinedRaw.raw"


awk 'NR==1 {
	for (i=1; i<=NF; i++) {
		sub(/_[^_]*$/, "", $i)
	}
} { print }' OFS='\t' "${PHENO_PATH}/testCombinedRaw.raw" > "${PHENO_PATH}/testCombined.raw"


rm "${PHENO_PATH}/testCombinedRaw.raw"


# Validate pheno_config file exists 
if [ ! -f "$PHENO_PATH/pheno.config" ]; then
	echo "Creating pheno_config file '${PHENO_PATH}/pheno.config' ... "
	touch "$PHENO_PATH/pheno.config"
	
else
	echo "Pheno_config file '${PHENO_PATH}/pheno.config' already exists ... "
fi


PHENO_CONFIG="$PHENO_PATH/pheno.config"

TEST_PATH="${PHENO_PATH}/testCombined.raw"
export TEST_PATH
# Replace line if TEST_PATH exists, else append it
if grep -q "^TEST_PATH=" "$PHENO_CONFIG"; then
	if [[ "$(uname)" == "Darwin" ]]; then
		# macOS sed syntax (requires '' for in-place)
		sed -i '' "s|^TEST_PATH=.*|TEST_PATH=${TEST_PATH}|" "$PHENO_CONFIG"
	else
		# Linux sed syntax
		sed -i "s|^TEST_PATH=.*|TEST_PATH=${TEST_PATH}|" "$PHENO_CONFIG"
	fi
else
	echo "TEST_PATH=${TEST_PATH}" >> "$PHENO_CONFIG"
	echo "export TEST_PATH" >> "$PHENO_CONFIG"

fi

TRAINING_PATH="${PHENO_PATH}/trainingCombined.raw"
export TRAINING_PATH
# Replace line if TRAINING_PATH exists, else append it
if grep -q "^TRAINING_PATH=" "$PHENO_CONFIG"; then
	if [[ "$(uname)" == "Darwin" ]]; then
		# macOS sed syntax (requires '' for in-place)
		sed -i '' "s|^TRAINING_PATH=.*|TRAINING_PATH=${TRAINING_PATH}|" "$PHENO_CONFIG"
	else
		# Linux sed syntax
		sed -i "s|^TRAINING_PATH=.*|TRAINING_PATH=${TRAINING_PATH}|" "$PHENO_CONFIG"
	fi
else
	echo "TRAINING_PATH=${TRAINING_PATH}" >> "$PHENO_CONFIG"
	echo "export TRAINING_PATH" >> "$PHENO_CONFIG"
fi

HOLDOUT_PATH="${PHENO_PATH}/holdoutCombined.raw"
export HOLDOUT_PATH
# Replace line if HOLDOUT_PATH exists, else append it
if grep -q "^HOLDOUT_PATH=" "$PHENO_CONFIG"; then
	if [[ "$(uname)" == "Darwin" ]]; then
		# macOS sed syntax (requires '' for in-place)
		sed -i '' "s|^HOLDOUT_PATH=.*|HOLDOUT_PATH=${HOLDOUT_PATH}|" "$PHENO_CONFIG"
	else
		# Linux sed syntax
		sed -i "s|^HOLDOUT_PATH=.*|HOLDOUT_PATH=${HOLDOUT_PATH}|" "$PHENO_CONFIG"
	fi
else
	echo "HOLDOUT_PATH=${HOLDOUT_PATH}" >> "$PHENO_CONFIG"
	echo "export HOLDOUT_PATH" >> "$PHENO_CONFIG"

fi

