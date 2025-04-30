#!/bin/bash
#
#SBATCH --job-name=create_filtered_snplist
#SBATCH -o /nfs/scratch/multerke/ukbiobank/err_out/%A_create_filtered_snplist.out
#SBATCH -e /nfs/scratch/multerke/ukbiobank/err_out/%A_create_filtered_snplist.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=2GB
#


cat "${PHENO_PATH}"/chr1_cleaned.snplist "${PHENO_PATH}"/chr2_cleaned.snplist "${PHENO_PATH}"/chr3_cleaned.snplist "${PHENO_PATH}"/chr4_cleaned.snplist \
"${PHENO_PATH}"/tanigawaSet/chr5_cleaned.snplist "${PHENO_PATH}"/tanigawaSet/chr6_cleaned.snplist "${PHENO_PATH}"/tanigawaSet/chr7_cleaned.snplist "${PHENO_PATH}"/chr8_cleaned.snplist \
"${PHENO_PATH}"/chr9_cleaned.snplist "${PHENO_PATH}"/chr10_cleaned.snplist "${PHENO_PATH}"/chr11_cleaned.snplist "${PHENO_PATH}"/chr12_cleaned.snplist \
"${PHENO_PATH}"/chr13_cleaned.snplist "${PHENO_PATH}"/chr14_cleaned.snplist "${PHENO_PATH}"/chr15_cleaned.snplist "${PHENO_PATH}"/chr16_cleaned.snplist \
"${PHENO_PATH}"/chr17_cleaned.snplist "${PHENO_PATH}"/chr18_cleaned.snplist "${PHENO_PATH}"/chr19_cleaned.snplist "${PHENO_PATH}"/chr20_cleaned.snplist \
"${PHENO_PATH}"/chr21_cleaned.snplist "${PHENO_PATH}"/chr22_cleaned.snplist "${PHENO_PATH}" > filteredSnps.snplist