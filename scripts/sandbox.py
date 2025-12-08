#!/usr/bin/env python3
import pandas as pd
import numpy as np


def find_opposite_sign_features(df, feature_sets,threshold):
    """
    Find features where the matching cohort has positive SHAP z-score 
    and ALL other cohorts have negative z-scores.
    
    Parameters:
    -----------
    df : DataFrame
        Must contain columns: 'features_used', 'cohort', 'feature', 'z_score', 'coef'
    feature_sets : list
        List of feature sets to analyze
    
    Returns:
    --------
    DataFrame with columns: feature_set, feature, matching_z_score, 
                           n_other_cohorts, coef, odds_ratio
    """
    rows = []
    
    for f in feature_sets:
        # Filter for current feature set
        df_feature_set = df[df['features_used'] == f]
        
        # Separate matching cohort from others
        df_matching = df_feature_set[(df_feature_set['cohort'] == f) & (df_feature_set['z_score'] > threshold)]
        df_other = df_feature_set[(df_feature_set['cohort'] != f) & (df_feature_set['feature'].isin(df_matching.feature.tolist()))]
        
        # Get unique features for this matching cohort
        for feature in df_matching['feature'].unique():
            # Get z-score and coef for this feature in matching cohort
            matching_row = df_matching[df_matching['feature'] == feature]
            
            if len(matching_row) == 0:
                continue
            
            matching_z = matching_row['z_score'].values[0]
            matching_coef = matching_row['coefs'].values[0]
            
            # Get z-scores for this feature in all other cohorts
            other_z = df_other[df_other['feature'] == feature]['z_score'].values
            
            # Skip if this feature doesn't appear in other cohorts
            if len(other_z) == 0:
                continue
            
            # Only keep if matching cohort has positive z-score 
            # AND ALL other cohorts have negative z-scores
            if matching_z > 0 and all(z < 0 for z in other_z):
                rows.append({
                    'feature_set': f,
                    'feature': feature,
                    'matching_z_score': matching_z,
                    'n_other_cohorts': len(other_z),
                    'coefs': matching_coef,
                    'odds_ratio': np.exp(matching_coef)
                })
                
    return pd.DataFrame(rows)

    
    
df = pd.read_csv('/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi/scores/importantFeaturesWithEpiMain/importantFeaturesAcrossCohortsAndTrainingDataThesisScript.csv')
#remove epi+main data
dfFilter1 = df[df['features_used'] != 'epi+main']
dfFiltered = dfFilter1[dfFilter1['cohort'] != 'epi+main']
training_data_list = ['HighCases','LowControls']

filteredResults = pd.DataFrame()

for training_data in training_data_list:
    df1 = dfFiltered[dfFiltered['training_data_used'] == training_data]
    feature_sets = df1['features_used'].unique()

    # Usage
    results_df = find_opposite_sign_features(df1, feature_sets,.5)
    results_df['training_data_used'] = training_data
    
    # Sort by odds ratio (descending) to see strongest effects first
    results_df = results_df.sort_values('odds_ratio', ascending=False)
    
    filteredResults = pd.concat([filteredResults,results_df],ignore_index=True)
    
filteredResults.to_csv('/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi/scores/importantFeaturesWithEpiMain/uniqueFeaturesAcrossCohortsAndTrainingData.Filtered.csv',index=False)
    


    