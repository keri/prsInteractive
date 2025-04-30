#!/bin/bash


plink --bfile "${PHENO_PATH}"/holdoutCombined --recode A --keep "${PHENO_PATH}/holdoutID.txt" --out "${PHENO_PATH}"/holdoutCombined


plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}/testID.txt" --out "${PHENO_PATH}"/testCombined


plink --bfile "${PHENO_PATH}"/merged_allChromosomes --recode A --keep "${PHENO_PATH}/trainingID.txt" --out "${PHENO_PATH}"/trainingCombined

#need one cleaned columns file to replace columns of .raw test and training set
python "${SCRIPTS_DIR}/clean_raw_columns01.py"

###############################################  Training Dataset ############################################

#remove the heading from the raw file into a temporary file
tail -n +2 "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/trainingCombined.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/trainingCombined_temp.raw

#take the columns file and insert as the first row of temp file and save as new file
cat "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/merged_allChromosomes_columns.txt "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/trainingCombined_temp.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/trainingCombined_final.raw

rm "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/trainingCombined_temp.raw

##############################################  Test Dataset ############################################

#remove the heading from the raw file into a temporary file
tail -n +2 "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/testCombined.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/testCombined_temp.raw

#take the columns file and insert as the first row of temp file and save as new file
cat "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/merged_allChromosomes_columns.txt "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/testCombined_temp.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/testCombined_final.raw

rm "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/testCombined_temp.raw

###############################################	 holdout Dataset	############################################ 

#remove the heading from the raw file into a temporary file
tail -n +2 "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/holdoutCombined.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/holdoutCombined_temp.raw

#take the columns file and insert as the first row of temp file and save as new file
cat "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/merged_allChromosomes_columns.txt "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/holdoutCombined_temp.raw > "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/holdoutCombined_final.raw

rm "${UKB_FILEPATH}"/"${pheno}"/tanigawaSet/data/holdoutCombined_temp.raw
