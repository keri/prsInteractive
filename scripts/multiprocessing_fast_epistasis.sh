#!/bin/bash



#Keri Multerer October 2022
#run fast fast-epistasis to generate seed snps


echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] number of cpus to run epistasis is set to: $N_CORES"


if [ ! -d "${PHENO_PATH}" ]; then
	echo "if you are running this script outside of the workflow you can either :"
	echo "2) COMMENT OUT PRECEDING STEPS IN RUN DATA_CLEANING_WORKFLOW OR"
	echo "3) run this as an isolated step using plink command"
	echo "for i in $(seq 1 ${N}); do
		cat output/directory/trainingEpi.epi.cc.${i} >> output/directory/trainingEpi.epi.c
		done
		plink --epistasis-summary-merge $OUTPUT_FILE $N --out output/directory/epiFiles/trainingCombinedEpi"
	echo "3) ./envSetUp.sh <pheno> <pheno string> <icd10> <n cores>"
	echo "   ./run_data_cleaning_workflow.sh <pheno> or on hpc: sbatch run_data_cleaning_workflow_submit.sh <pheno> "
	
else
	#EPI_PATH is exported from data cleaning workflow
	OUT_DIR="${EPI_PATH}/preSummaryFiles"
	OUTPUT_FILE="${OUT_DIR}/trainingEpi.epi.cc"
fi 


# Loop through parallel chunks 1 to 40
for i in $(seq 1 ${N_CORES});
do
	echo "in the plink loop : "${i}
	plink --bfile "${PHENO_PATH}"/merged_allChromosomes --keep "${PHENO_PATH}/trainingID.txt" --fast-epistasis 'boost' --parallel $i $N_CORES --out "${OUT_DIR}"/trainingEpi
	
done

# Wait for all background jobs to finish
wait

echo "All $N_CORES parallel PLINK jobs finished."

# Empty the output file if it already exists
> $OUTPUT_FILE

for i in $(seq 1 ${N_CORES}); 
do
	cat "${OUT_DIR}/trainingEpi.epi.cc.${i}" >> $OUTPUT_FILE
done

plink --epistasis-summary-merge $OUTPUT_FILE $N_CORES --out "${PHENO_PATH}/epiFiles/trainingCombinedEpi"

PHENO_CONFIG="$PHENO_PATH/pheno.config"

# Validate pheno_config file exists 
if [ ! -f "$PHENO_CONFIG" ]; then
	echo "Creating pheno_config file '$PHENO_CONFIG' ... "
	touch "$PHENO_CONFIG"
	
else
	echo "Pheno_config file '$PHENO_CONFIG' already exists ... "
fi

EPI_FILE="$PHENO_PATH/epiFiles/trainingCombinedEpi.epi.cc.summary"
echo "Epi path is set to ... $EPI_FILE"

#check to see if EPI_PATH exists
if grep -q "^EPI_FILE=" "$PHENO_CONFIG"; then
	if [[ "$(uname)" == "Darwin" ]]; then
		# macOS sed syntax (requires '' for in-place)
		sed -i '' "s|^EPI_FILE=.*|EPI_FILE=${EPI_FILE}|" "$PHENO_CONFIG"
	else
		# Linux sed syntax
		sed -i "s|^EPI_FILE=.*|EPI_FILE=${EPI_FILE}|" "$PHENO_CONFIG"
	fi
else
	echo "EPI_FILE=${EPI_FILE}" >> "$PHENO_CONFIG"
	echo "export EPI_FILE" >> "$PHENO_CONFIG"	
fi

