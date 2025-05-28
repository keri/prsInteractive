#!/bin/bash
#Keri Multerer October 2022
#run fast fast-epistasis to generate seed snps


echo "[BASH] Reading from: $PHENO_PATH"
echo "[WORKFLOW] number of cpus to run epistasis is set to: $N"

OUT_DIR="${PHENO_PATH}/epiFiles/preSummaryFiles"
OUTPUT_FILE="${OUT_DIR}/trainingEpi.epi.cc"

if [ ! -d "${OUT_DIR}" ]; then
	echo "Folder '${PHENO_PATH}'/epiFiles/preSummaryFiles does not exist. Creating it..."
	mkdir "${PHENO_PATH}/epiFiles" 
	mkdir "${PHENO_PATH}/epiFiles/preSummaryFiles" 
	
else
	echo "Folder $OUT_DIR already exists."	
fi


# Loop through parallel chunks 1 to 40
for i in $(seq 1 ${N});
do
	echo "in the plink loop : "${i}
	plink --bfile "${PHENO_PATH}"/merged_allChromosomes --keep "${PHENO_PATH}/trainingID.txt" --fast-epistasis 'boost' --parallel $i $N --out "${OUT_DIR}"/trainingEpi
	
done

# Wait for all background jobs to finish
wait

echo "All 40 parallel PLINK jobs finished."

# Empty the output file if it already exists
> $OUTPUT_FILE

for i in $(seq 1 ${N}); 
do
	cat "${OUT_DIR}/trainingEpi.epi.cc.${i}" >> $OUTPUT_FILE
done

plink --epistasis-summary-merge $OUTPUT_FILE $N --out "${PHENO_PATH}/epiFiles/trainingCombinedEpi"

PHENO_CONFIG="$PHENO_PATH/pheno_config.sh"

# Validate pheno_config file exists 
if [ ! -f "$PHENO_CONFIG" ]; then
	echo "Creating pheno_config file '$PHENO_CONFIG' ... "
	touch "$PHENO_CONFIG"
	
else
	echo "Pheno_config file '$PHENO_CONFIG' already exists ... "
fi

EPI_PATH="$PHENO_PATH/epiFiles/trainingCombinedEpi.epi.cc.summary"
export EPI_PATH
#check to see if EPI_PATH is already present
if grep -q "^EPI_PATH=" "$PHENO_CONFIG"; then
	if [[ "$(uname)" == "Darwin" ]]; then
		# macOS sed syntax (requires '' for in-place)
		sed -i '' "s|^EPI_PATH=.*|EPI_PATH=${EPI_PATH}|" "$PHENO_CONFIG"
	else
		# Linux sed syntax
		sed -i "s|^EPI_PATH=.*|EPI_PATH=${EPI_PATH}|" "$PHENO_CONFIG"
	fi
else
	echo "EPI_PATH=$EPI_PATH" >> "$PHENO_CONFIG"
	echo "export EPI_PATH" >> "$PHENO_CONFIG"
	echo "echo EPI_PATH is set to: $EPI_PATH" >> "$PHENO_CONFIG"

fi



	