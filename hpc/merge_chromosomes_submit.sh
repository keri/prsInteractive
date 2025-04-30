#!/bin/bash
#
#SBATCH --job-name=merge_chromosomesTanigawa
#SBATCH -o  /nfs/scratch/multerke/ukbiobank/err_out/%A_mergeTrainHoldout.out
#SBATCH -e /nfs/scratch/multerke/ukbiobank/err_out/%A_mergeTrainHoldout.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=4:00:00
#

module load plink/1.90


plink --bfile "${DATA_PATH}"/chr1_cleaned --merge-list "${DATA_PATH}"/mergeChromosomes.txt  --write-snplist --make-bed --out "${DATA_PATH}"/merged_allChromosomes

plink --bfile "${DATA_PATH}"/variant_calls/ukb22418_c1_b0_v2 --merge-list "${DATA_PATH}"/variant_calls/mergeChomosomesHoldout.txt --keep "${PHENO_PATH}"/holdoutID.txt --pheno "${PHENO_PATH}"/pheno.txt --extract "${DATA_PATH}"/merged_allChromosomes.snplist --make-bed --out "${PHENO_PATH}"/holdoutCombined
