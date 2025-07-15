import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import warnings
warnings.simplefilter(action='ignore')
import sys
import os
# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
from .download import *
from .draw_plots import *

def create_box_plot_covar(prs,figurePath,filepath):
    prsCopy = prs.copy()
    data_type = 'covariate'
    prsCopy['prs'] = prsCopy.drop(columns=['PHENOTYPE']).sum(axis=1)

    scaler = StandardScaler().set_output(transform="pandas")
    scaled_feature = scaler.fit_transform(prsCopy[['prs']])


    prsCopy['scaled_prs'] = scaled_feature

    prsCopy[['PHENOTYPE','prs','scaled_prs']].to_csv(f'{filepath}prs/prs_{data_type}.csv')

    ax = prsCopy.boxplot(column='prs',by='PHENOTYPE',figsize=(10,10))
    ax.plot()
    title = f'prs for dataset {data_type}'
    plt.title(title)
    plt.savefig(f'{figurePath}{data_type}_prs.png')
    plt.close()


def create_prs_plots(prs,data_type,figurePath,filePath,direction):
    prsCopy = prs.copy()

    str_image = f'{data_type}.{direction}'
    prs.to_csv(f'{filePath}/{str_image}.prs.csv')
    create_box_plot(prsCopy, str_image, figurePath)
    create_density_plot(prsCopy,str_image,figurePath)
    threshold = create_auc_graph(prsCopy,str_image,figurePath)
    create_prevalence_plot(prsCopy,str_image,figurePath)
    
def get_feature_coef_dictionaries(featureScores):
    featureScoresTemp = featureScores.copy()
    featureScoresTemp.set_index(['feature'],inplace=True)
    featureScoresDict = featureScoresTemp[['coefs']].T.to_dict()
    return(featureScoresDict)

def standardize_prs(prsDf):
    prsCopy = prsDf.copy()
    prsCopy['prs'] = prsCopy.drop(columns=['PHENOTYPE']).sum(axis=1)
    
    scaler = StandardScaler().set_output(transform="pandas")
    scaled_feature = scaler.fit_transform(prsCopy[['prs']])
    
    prsCopy['scaled_prs'] = scaled_feature
    return(prsCopy)

def calculate_prs(df,featureScoreDict,features):
    '''calculate PRS for every person using beta coefficients from regression stored in dictionary
       input : nested dictionary {feature1: {coefs:beta},.... featureN:{coefs:beta}}
               features : list of strings '''     
    dfCopy = df.copy()
    prsDf = dfCopy[['PHENOTYPE']]

    #if feature is an epi pair, add columns and create feature in prsCombined, otherwise ignore and add as is
    for feature in features:
#       dfFeature = pd.DataFrame()
#       col = feature.split(',')
#       dfFeature[feature] = dfCopy[col].sum(axis=1)
        try:
#           prsDf[feature] = dfFeature[feature].apply(lambda x : (featureScoreDict[feature]['coefs'])*x)
            prsDf[feature] = dfCopy[feature].apply(lambda x : (featureScoreDict[feature]['coefs'])*x)
        except KeyError:
            print(f'{feature} not in data set ..')
            
    prsDf2 = standardize_prs(prsDf)
    try:
        meanDiff = prsDf2.groupby(['PHENOTYPE']).mean()['scaled_prs'].diff().loc[2]
        print('mean diff after box plot created ',meanDiff)
        
    except KeyError:
        meanDiff = np.nan
        print('unable to calculate mean diff of the 2 populations, possibly due to insufficient numbers in 1 population (cases)  ',meanDiff)
        
    return(prsDf2,meanDiff)

def calculate_create_prs_plots(df,featureScoreDict,data_type,figurePath,prsPath,direction,features):
    '''df dataframe : columns['PHENOTYPE',str(var1),str(var2)...str(varN)
    featureScoreDict : dict key=variant, value=coef
    data_type : str(main,epi,epi+main)
    figurePath : str(folder path to save figures)
    prsPath : str(folder to save calculated prs scores columns[IID,PHENOTYPE,prs,scaled_prs]
    direction : str(risk,protect,or mixed)
    features : list(features to include in prs calculation)'''
    
    prsDf,meanDiff = calculate_prs(df,featureScoreDict,features)
#   if 'holdout' in data_type:
#       #switch the prs and scaled_prs columns to plot prs instead of scaled prs without refactoring all functions
#       prsDf.rename(columns={'scaled_prs':'scaled_prs_temp'},inplace=True)
#       prsDf.rename(columns={'prs':'scaled_prs'},inplace=True)
#       prsDf.rename(columns={'scaled_prs_temp':'prs'},inplace=True)
        
    create_prs_plots(prsDf,data_type,figurePath,prsPath,direction)
    return(meanDiff)
    
        
    



    