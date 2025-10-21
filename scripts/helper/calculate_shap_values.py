#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
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
    
    if data_type == 'topFeatures': #this is for the top feature analysis (27_calculate_top_features_in_PRSiPRS_group.py)
#       topFeatures = stats.zscore(shap_values3)
        topFeatures = shap_values4[(shap_values4 > 2) | (shap_values4 < -2) ]
#       topFeatures = shap_values4[shap_values4 > 1].sort_values(ascending=True)
    elif data_type == 'topFeatures.HighCases':
        topFeatures = shap_values4[(shap_values4 > 2) | (shap_values4 < -2) ]
    else:#this can be changed for some of the epi features which are much higher than main for some phenotypes
        topFeatures = shap_values4[(shap_values4 > 2) | (shap_values4 < -2)].sort_values(ascending=True)
    
    topFeaturesMean = shap_values3.loc[topFeatures.index]
    return(topFeatures,shap_values3,shap_values4)

def create_summary_plots(figPath,shap_values,X_test,i,data_type):
    # Create the SHAP summary plot
    plt.figure(figsize=(15,10))  # Start a new figure
    shap.summary_plot(shap_values, X_test, show=False)
    
    # Save the plot as a file
    plt.savefig(f"{figPath}/shap_summary_plot.{i}.{data_type}.png", dpi=300, bbox_inches="tight")
    plt.close()  # Close the figure
    
    
    
def plot_and_save_top_features(featuresZscore,i,figPath,data_type):
#   features = topFeatures.index.tolist()
    plt.figure(figsize=(8, 10))
    try:
        featuresZscore.plot(kind='barh',figsize=(10,20),title=f'z scores using {data_type} for {i}')
        plt.savefig(f'{figPath}/importantFeatureZscores.{i}.{data_type}.png',format="png", dpi = 300, bbox_inches='tight')
        plt.close()
    except IndexError:
        pass
        
#def calculate_plot_shap_values(model,X_test,y_test,i,figPath,data_type):
#   
#   #features = X_test.columns.tolist()
#   clfHGB = model.best_estimator_
#   features = clfHGB.feature_names_in_
#   explainer = fasttreeshap.TreeExplainer(clfHGB, algorithm='auto',n_jobs=-1)
#   shap_explainer = explainer(X_test,check_additivity=False)
#   shap_values = shap_explainer.values
#   index = X_test.index
#   
#   #shap_file = f"{trainingPath}/featureShapValue_{data_type}_{i}.csv"
#   #get the shap values in dataframe form
#   shap_df = pd.DataFrame(data=shap_values, columns=features, index=index)
#   #shap_df.to_csv(shap_file)
#   
#   topFeatures,featureMean,shap_valuesZscores = get_top_shap_features(shap_df,data_type)
#   
#   create_summary_plots(figPath,shap_values,X_test,i,data_type)
#   
#   plot_and_save_top_features(topFeatures,i,figPath,data_type)
#   
#   return(topFeatures,shap_valuesZscores)

def calculate_plot_shap_values(model, X_test, y_test, i, figPath, data_type):
    
    # Get the best XGBoost estimator
    clfXGB = model.best_estimator_
    
    # Check if the attribute exists before accessing it
    if hasattr(model, 'feature_names_in_'):
        feature_names = model.feature_names_in_
    else:
        # Fallback: use column names from your DataFrame or create generic names
        feature_names = X_test.columns.tolist()  # if X_train is a DataFrame
        # OR create generic names:
        # feature_names = [f'feature_{i}' for i in range(X_train.shape[1])]
        
#   print(f"Feature names: {feature_names[:10]}")

    # Use shap with XGBoost model

        #explainer = shap.TreeExplainer(clfXGB)
    try:
        explainer = shap.TreeExplainer(
            clfXGB, 
            feature_perturbation='interventional',
            model_output='raw'
        )
        shap_values = explainer.shap_values(X_test)
        print("✓ Standard SHAP with interventional worked")
    except Exception as e:
        print(f"✗ Standard SHAP failed: {e}")

    
        # Option 2: Use Permutation explainer as fallback
        try:
            explainer = shap.PermutationExplainer(clfXGB.predict, X_test[:100])
            shap_values = explainer.shap_values(X_test)
            print("✓ Permutation explainer worked (slower but reliable)")
    
        except Exception as e:
            print(f"✗ Permutation explainer failed: {e}")
            return None, None

#   explainer = fasttreeshap.TreeExplainer(clfXGB)
#   shap_values = explainer.shap_values(X_test)
    
    # Handle multi-class output (shap_values might be a list)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]  # For binary classification, use positive class
    
    index = X_test.index
    
    # Get the shap values in dataframe form
    shap_df = pd.DataFrame(data=shap_values, columns=feature_names, index=index)
    
    topFeatures, featureMean, shap_valuesZscores = get_top_shap_features(shap_df, data_type)
    
    create_summary_plots(figPath, shap_values, X_test, i, data_type)
    
    plot_and_save_top_features(topFeatures, i, figPath, data_type)
    
    return(topFeatures, shap_valuesZscores)





    