#!/usr/bin/env python3

import pandas as pd
from sklearn.preprocessing import StandardScaler
import os
import argparse

def rank_gene_env_features(geneEnvShapleyFile,threshold=2):
    '''
    input : df with columns [envGeneticFeature,shap_zscore,env_type,geneticFeature,envFeature,main_E,epistatic]
    
    oupt : df with gene_environment_features ranked with Shapley Z scores > 2
        
    '''
    #download data
    df = pd.read_csv(geneEnvShapleyFile)
    
    #get the largest shap_zscore for duplicated values
    df.sort_values(['envGeneticFeature','shap_zscore'],ascending=False,inplace=True)
    
    #drop duplicates in df
    df.drop_duplicates(subset=['envGeneticFeature'],keep='first',inplace=True)
    
    #get the epistatic interactions
    epiDf = df[df['epistatic'] == 1]
    
    #get the main effects
    mainDf = df[df['main_E'] == 1]
    mainDf.loc[mainDf['envFeature'].isna(),'envFeature'] = mainDf['envGeneticFeature']
    
    importantFeatures = epiDf[epiDf['shap_zscore'] > 3]
    
    finalFeatures = pd.concat([importantFeatures,mainDf],ignore_index=True)
    
    newFile = geneEnvShapleyFile.split('.')[0]
    newFile = f'{newFile}Full.csv'
    df.to_csv(newFile,index=False)
    
    finalFeatures.to_csv(geneEnvShapleyFile,index=False)
    
    return finalFeatures

def main(featureFile):
    
    df = pd.read_csv(featureFile)
    
    reducedFeatures = rank_gene_env_features(df)
    
    newFile = featureFile.split('.')[0]
    newFile = f'{newFile}Full.csv'
    df.to_csv(newFile,index=False)
    
    reducedFeatures.to_csv(featureFile,index=False)

if __name__ == '__main__':
    

    filePath = '/Users/kerimulterer/prsInteractive/results/myocardialInfarction/scores/cardioMetabolicimportantFeaturesPostShap.csv'

    importantFeatures = rank_gene_env_features(filePath)
        