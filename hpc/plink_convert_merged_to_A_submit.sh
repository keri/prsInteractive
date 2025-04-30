#!/bin/bash
#
#SBATCH --job-name=convert_merged_dataset_to_raw
#SBATCH -o  /nfs/scratch/multerke/ukbiobank/err_out/%A_convert_merged_dataset_to_raw.out
#SBATCH -e /nfs/scratch/multerke/ukbiobank/err_out/%A_convert_merged_dataset_to_raw.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=60G
#SBATCH --time=05:00:00
#


module load plink/1.90

plink --bfile "${PHENO_PATH}"/holdoutCombined --recode A --keep "${PHENO_PATH}/holdoutID.txt" --out "${PHENO_PATH}"/holdoutCombined


plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}/testID.txt" --out "${PHENO_PATH}"/testCombined


plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}/trainingID.txt" --out "${PHENO_PATH}"/trainingCombined

module load Miniconda3/4.9.2
eval "$(conda shell.bash hook)"
conda activate /nfs/home/multerke/.conda/envs/py39

#need one cleaned columns file to replace columns of .raw test and training set
sbatch 09A_clean_raw_columns_submit.sh

#replaces the column file for both test and training data
sbatch 09B_use_cleaned_raw_columns_submit.sh