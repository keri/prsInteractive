#!/bin/bash



#Keri Multerer October 2022
#run fast fast-epistasis to generate seed snps


echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"

N=18

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

	