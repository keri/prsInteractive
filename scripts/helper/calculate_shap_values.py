#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shap
import fasttreeshap
shap.initjs()




def get_top_shap_features(shap_values,data_type):
    '''standardize shap values and get features with Z score > 2 
    returns list of shap features with z scores > threshold'''
    
    #remove the features that have value of 0
    shap_values2 = shap_values.loc[:,(shap_values != 0).any(axis=0)]
    
    #get the mean for each remaining feature
#   shap_values3 = shap_values2.apply(np.abs).mean().sort_values(ascending=False)
    shap_values3 = shap_values2.mean().sort_values(ascending=False)
    
    #get the z score for the feature means
    shap_values4 = stats.zscore(shap_values3)
    if data_type == 'main':
        topFeatures = shap_values4[shap_values4 > 2].sort_values(ascending=True)
    elif data_type == 'topFeatures': #this is for the top feature analysis (27_calculate_top_features_in_PRSiPRS_group.py)
#       topFeatures = stats.zscore(shap_values3)
        topFeatures = shap_values4[(shap_values4 > 2) | (shap_values4 < -2) ]
#       topFeatures = shap_values4[shap_values4 > 1].sort_values(ascending=True)
    elif data_type == 'topFeatures.HighCases':
        topFeatures = shap_values4[(shap_values4 > 2) | (shap_values4 < -2) ]
    elif data_type == 'cardioMetabolic':
        topFeatures = shap_values4[shap_values4 > 2].sort_values(ascending=True)
    elif data_type == 'mentalHealthNutrition':
        topFeatures = shap_values4[shap_values4 > 2].sort_values(ascending=True)
    else:#this can be changed for some of the epi features which are much higher than main for some phenotypes
        topFeatures = shap_values4[shap_values4 > 5].sort_values(ascending=True)
    
    
    return(topFeatures,shap_values3,shap_values4)

def create_summary_plots(figPath,shap_values,X_test,i):
    # Create the SHAP summary plot
    plt.figure(figsize=(15,10))  # Start a new figure
    shap.summary_plot(shap_values, X_test, show=False)
    
    # Save the plot as a file
    plt.savefig(f"{figPath}/shap_summary_plot.{i}.png", dpi=300, bbox_inches="tight")
    plt.close()  # Close the figure



def plot_and_save_top_features(featuresZscore,i,figPath,data_type):
#   features = topFeatures.index.tolist()
    try:
        featuresZscore(kind='barh',figsize=(10,20),title=f'z scores using {data_type} for {i}')
        plt.savefig(f'{figPath}/featureZScores.{i}.{data_type}.png',format="png", dpi = 300, bbox_inches='tight')
        plt.close()
    except IndexError:
        pass

def calculate_plot_shap_values(model,X_test,y_test,i,figPath,data_type):
    
    features = X_test.columns.tolist()
    clfHGB = model.best_estimator_
    explainer = fasttreeshap.TreeExplainer(clfHGB, algorithm='auto',n_jobs=-1)
    shap_explainer = explainer(X_test,check_additivity=False)
    shap_values = shap_explainer.values
    index = df.index
    
    #shap_file = f"{trainingPath}/featureShapValue_{data_type}_{i}.csv"
    #get the shap values in dataframe form
    shap_df = pd.DataFrame(data=shap_values, columns=features, index=index)
    #shap_df.to_csv(shap_file)
    
    shap_features,topFeatures,shap_valuesMean,shap_valuesZscores = get_top_shap_features(shap_df,data_type)
    
    create_summary_plots(figPath,shap_values,X_test,i)
    
    plot_and_save_top_features(shap_valuesZscores,i,figPath,data_type)
    
    return(topFeatures)

    