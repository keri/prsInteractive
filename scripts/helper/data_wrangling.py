#!/usr/bin/env python3

import pandas as pd
import numpy as np

def create_epi_df(epiDf,pairList):
    '''input : epiDf with snps as columns + PHENOTYPE
            pairList : [snp1pair1,snp2pair1,snp1pair2,snp2pair2...snp1pairN,snp2pairN]'''
    epiArrayFinal = pd.DataFrame()
    
    for pair in pairList:
        snps = epiDf[pair.split(',')].sum(axis=1)
        snps.columns = pair
        epiArrayFinal = pd.concat([epiArrayFinal,snps],axis=1)
    epiArrayFinal.columns = pairList 
    
    return(epiArrayFinal)