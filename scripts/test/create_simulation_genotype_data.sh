#!/bin/bash

for ((a=1; a<23; a++)); do
    plink --simulate "${SCRIPTS_DIR}"/test/gwas.sim --simulate-ncases 500 --simulate-ncontrols 1500 --make-bed --out "${DATA_PATH}"/variant_calls/ukb22418_c"${a}"_b0_v2
done
