#!/bin/bash

for ((a=1; a<23; a++)); do
    plink --simulate "${SCRIPTS_DIR}"/test/gwas.sim --simulate-ncases 500 --simulate-ncontrols 1500 --make-bed --out "${DATA_PATH}"/variant_calls/temp
    
    # Remove first row AND "per" prefix
    tail -n +2 "${DATA_PATH}"/variant_calls/temp.fam | sed 's/per//g' > "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2.fam
    
    # Copy other files
    cp "${DATA_PATH}"/variant_calls/temp.bim "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2
    cp "${DATA_PATH}"/variant_calls/temp.bed "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2
    
    # Clean up
    rm "${DATA_PATH}"/variant_calls/temp.*
done
