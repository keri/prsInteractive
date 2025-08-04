#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

pheno = 'type2Diabetes'
scoresPath = f'/Users/kerimulterer/prsInteractive/results/{pheno}/scores'
featureFile = f'{scoresPath}/featureScoresReducedFinalModel.filtered.csv'
hlaFile = '/Users/kerimulterer/prsInteractive/results/participant_hla.csv'

def get_hla_columns(hlaFile):
    df = pd.read_csv(hlaFile,nrows=1)
    try:
        df.drop(columns=['IID'],inplace=True)
        hlaColumns = df.columns.tolist()
    except KeyError:
        hlaColumns = df.columns.tolist()
        
    return hlaColumns

def calculate_mean_for_features(df):
    return df[['coefs']].mean(axis=0)

def graph_mean_diffs(meanDiffsDict):
    

def calculate_mean_diff_across_models(scoresPath,hlaFile,featureFile):
    
    featureDf = pd.read_csv(featureFile)
    
    covarFeatures = featureDf[featureDf['model'] == 'covariate']['feature'].tolist()
    
    #get only relevant models
    featureDf = featureDf[featureDf['model'].isin(['main','epi','epi+main','cardio','all'])]
    
    #remove covariate feeatures
    featureDf = featureDf[~featureDf['feature'].isin(covarFeatures)]
    
    meanDiffs = {}
    #get the feature sets to compare with all feature
    for model in featureDf[featureDf['model'] != 'all']['model'].unique():
        feature_set = featureDf[featureDf['model'] == model]['feature'].tolist()
        all_feature_df = featureDf[featureDf['model'] == 'all']
        mean_all = calculate_mean_for_features(all_feature_df[all_feature_df['feature'].isin(feature_set)])
        mean_model = calculate_mean_for_features(featureDf[featureDf['model'] == model])
        meanDiffs[model] = mean_all-mean_model
        
    



