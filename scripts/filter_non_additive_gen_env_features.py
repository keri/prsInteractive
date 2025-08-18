#!/usr/bin/env python3

import pandas as pd
import os
import argparse

def compare_gene_env_to_genetic(scoresPath):
    inputFile = f'{scoresPath}/featureScoresReducedFinalModel.csv'
    df = pd.read_csv(inputFile)
    
    #columns = coefs,model,feature
    #iterate over all non-zero GxGxE features and get the first element in list for interactive features to compare with epi or main models
    gene_env = df[(df['model'] == 'cardio') & (df['coefs'] != 0) & (df['feature'].str.contains(','))]
    filterList = []
    for feature in gene_env['feature']:
        g = ','.join(feature.split(',')[1:])
        eDf = gene_env[gene_env['feature'] == feature]
        #get coef value for single g features
        gDf = df[(df['feature'] == g) & (df['model'].isin(['epi','main']))]
        if gDf.empty:
            pass
        else:
            gValue = gDf.sort_values(['coefs'],ascending=False).head(1)['coefs'].values[0]
            eValue = eDf['coefs'].values[0]
            if gValue > 0: #for GxGxE features that are more predictive of risk
                if gValue > eValue:
                    filterList.append(feature)
            else:
                if eValue > gValue: #for GxGxE features that more predictive of controls
                    filterList.append(feature)
                
    dfFiltered = df[~df['feature'].isin(filterList)]
    inputFilePrefix = inputFile.split('.')[0]
    dfFiltered.to_csv(f'{inputFilePrefix}.filtered.csv')
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="calculating PRS for association modelling weights....")
    parser.add_argument("--scores_path", help="Path to the input scores folder")
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    scores_path = args.scores_path or os.environ.get['SCORES_PATH']
        
    main(scores_path)