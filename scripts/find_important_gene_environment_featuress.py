#!/usr/bin/env python3

import pandas as pd
from sklearn.preprocessing import StandardScaler
import os
import argparse

def rank_gene_env_features(geneEnvShapleyFile,threshold):
    '''
    input : df with columns [envGeneticFeature,shap_zscore,env_type,geneticFeature,envFeature,main_E,epistatic]
    
    oupt : df with gene_environment_features ranked with Shapley Z scores > 2
        
    '''
    #download data
    df = pd.read_csv(geneEnvShapleyFile)
    
    #remove the main_E = 2 as that feature is not being used
    df1 = df[df['main_E'].isin([0,1])]
    
    #take the main of all main E features
    epiDf = df1[df1['epistatic'] == 1].groupby(by=['envGeneticFeature','env_type', 'geneticFeature',
        'envFeature', 'main_E', 'epistatic']).mean().reset_index()
    
    #get the epistatic interactions
#   epiDf = df2[df2['epistatic'] == 1]
    
    #get the main effects
    mainDf = df1[df1['main_E'] == 1]
    mainDf.loc[mainDf['envFeature'].isna(),'envFeature'] = mainDf['envGeneticFeature']
    
    mainDfCleaned = mainDf.groupby(by=['envGeneticFeature','env_type',
        'envFeature', 'main_E', 'epistatic']).mean().reset_index()
    
    
    
    
    #get the rows with important shap values
    importantFeatures = epiDf[(epiDf['shap_zscore'] > threshold) | (epiDf['shap_zscore'] < -threshold)]
    
#   #get the largest shap_zscore for duplicated values
#   importantFeatures.sort_values(['envGeneticFeature','shap_zscore'],ascending=False,inplace=True)
#   
#   #drop duplicates in df
#   importantFeatures.drop_duplicates(subset=['envGeneticFeature'],keep='first',inplace=True)
    
    
    finalFeatures = pd.concat([importantFeatures,mainDfCleaned],ignore_index=True)
    
    newFile = geneEnvShapleyFile.split('.')[0]
    newFile = f'{newFile}FilteredZscore.csv'
#   df.to_csv(newFile,index=False)
    
    finalFeatures.to_csv(newFile,index=False)
    
    return finalFeatures


def main(featureFile,threshold):
    
    df = pd.read_csv(featureFile)
    
    reducedFeatures = rank_gene_env_features(df,threshold)
    
#   newFile = featureFile.split('.')[0]
#   newFile = f'{newFile}Full.csv'
#   df.to_csv(newFile,index=False)
    
#   reducedFeatures.to_csv(featureFile,index=False)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="filtering important GxGxE features to use...")
    parser.add_argument("--gene_env_file",help="Genetic environmental epistatic features")
    parser.add_argument("--threshold", type=float, default=1.99, help="shapley z score threshold to use for filter (default: %(default)s)")
    
    args = parser.parse_args()
    
    gene_env_file = args.gene_env_file or os.environ.get("GENE_ENV_FILE")
    print(f"reading from gene environment epi features file : {gene_env_file}")
    
    threshold = args.threshold or os.environ.get("THRESHOLD")
    threshold = float(threshold)
    print(f"analyzing top features based on shap z score : {threshold}")
    
#   gene_env_file = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/productEpi/scores/cardioMetabolicimportantFeaturesPostShap.csv'
#   threshold = 1.99

    importantFeatures = rank_gene_env_features(gene_env_file,threshold)
        