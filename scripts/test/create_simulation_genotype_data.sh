#!/bin/bash

for ((a=1; a<23; a++)); 
do
    plink --simulate ~/prsInteractive/testData/wgas.sim --make-bed --out ~/prsInteractive/testData/sim_variant_calls/ukb22418_c"${a}"_b0_v2
    
done
