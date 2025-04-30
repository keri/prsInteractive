#!/bin/bash
#Keri Multerer October 4022
#run fast fast-epistasis to generate seed snps
#
#SBATCH --job-name=fast_fast-epistasis_tanigawa
#SBATCH -o  /nfs/scratch/multerke/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/multerke/ukbiobank/err_out/%A.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --time=62:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keri@multerer.com
#
module load plink/1.90

echo "[BASH] Reading from: $PHENO_PATH"
cat "$PHENO_PATH"

OUT_DIR="$PHENO_PATH/epiFiles/preSummaryFiles"
OUTPUT_FILE="${PHENO_PATH}/epiFiles/trainingCombinedEpi.epi.cc"

if [ ! -d "${OUT_DIR}/epiFiles" ]; then
	echo "Folder '${PHENO_PATH}/epiFiles/preSummaryFiles' does not exist. Creating it..."
	mkdir "${PHENO_PATH}/epiFiles" 
	mkdir "${PHENO_PATH}/epiFiles/preSummaryFiles" 
	
else
	echo "Folder '${PHENO_PATH}/epiFiles/preSummaryFiles' already exists."	
fi

# Loop through parallel chunks 1 to 40
for i in $(seq 1 40);
do
plink --bfile "${PHENO_PATH}/merged_allChromosomes" --keep "${PHENO_PATH}/trainingID.txt" --remove  "${DATA_PATH}/withdrawals.csv" --fast-epistasis 'boost' --parallel 1 ${i} --out "${OUT_DIR}/trainingEpi"

done

# Wait for all background jobs to finish
wait

echo "All 40 parallel PLINK jobs finished."

# Empty the output file if it already exists
> $OUTPUT_FILE

for i in $(seq 1 40); 
do
	cat "${OUT_DIR}/trainingEpi.epi.cc.${i}" >> $OUTPUT_FILE
done

plink --epistasis-summary-merge $OUTPUT_FILE 40 --out "${PHENO_PATH}/epiFiles/trainingCombinedEpi"

	