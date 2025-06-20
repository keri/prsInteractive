#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

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

def get_epi_snps(epiFeatures):
    '''input : list[str:pair1,str:pair2..str:pairN]
    output : list : unique[str:snp1,str:snp2...str:snpN]'''
    epiSnps = []
    for pair in epiFeatures:
        epiSnps = epiSnps + pair.split(',')
    epiSnps = list(set(epiSnps))
    return(epiSnps)

def scale_data(df):
    # Initialize the StandardScaler
    scaler = StandardScaler()
    
    # Fit the scaler to your data (compute the mean and standard deviation)
    scaler.fit(df)
    
    # Transform the data using the fitted scaler
    scaled_data = scaler.transform(df)
    
    # Create a new DataFrame with the scaled data
    scaled_df = pd.DataFrame(scaled_data, columns=df.columns,index=df.index)
    
    return(scaled_df)





#   
#   
#
#
#if __name__ == '__main__':
#   
#   main('/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/epiFiles')


    