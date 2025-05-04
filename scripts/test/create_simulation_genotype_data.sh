#!/bin/bash

for ((a=1; a<23; a++)); 
do
    plink --simulate ~/prsInteractive/testData/gwas.sim --make-bed --out ~/prsInteractive/testData/variant_calls/ukb22418_c"${a}"_b0_v2
    
done
