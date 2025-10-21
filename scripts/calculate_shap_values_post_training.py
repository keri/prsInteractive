#!/usr/bin/env python3

import pandas as pd
import numpy as np
import time
#from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.impute import SimpleImputer
import pickle
import csv
import os
import argparse
import glob

from helper.download import get_dataset, get_epi_columns, get_columns, get_column_index
from helper.calculate_shap_values import *
from helper.data_wrangling import *



def calculate_shap_features(X,y,pheno,data_type,modelFile,imp_mean,clfHGB,figPath): 
    '''load pickled models and score with test set'''
    print('scoring models .....')
    

def download_data(pheno_path,i,data_type,test_path,epi_path,withdrawal_path):
    full_columns = get_columns(test_path)
    n=3000

    #############################################################################
    #                 GET SNPS/SNP PAIRS FOR MODEL                              #
    ############################################################################# 
    

    istart = (i-1)*n
    if data_type == 'main':
#           if i == 0:
#               sectionSnps = full_columns[:istart+n]
#           else:
        
        sectionSnps = full_columns[istart+6:istart+n+6]
        sectionPairs = sectionSnps
        sectionPairs = [x for x in sectionPairs if pd.notna(x)]
        
    else:
        print('epi_path to get filtered epi pairs ..',epi_path)
        epiColumns = get_epi_columns(epi_path)
        sectionPairs = epiColumns[istart:istart+n]
        sectionPairs = [x for x in sectionPairs if pd.notna(x)]
        sectionSnps = get_epi_snps(sectionPairs)
    mainArray = get_dataset(test_path,withdrawal_path,sectionSnps, use_chunking=True)
    return(mainArray,sectionPairs)
            
    
    
def main(pheno_path,test_path,epi_file,withdrawal_path):
    figPath = f'{pheno_path}/figures'
    scoresPath = f'{pheno_path}/scores'
    modelPath = f'{pheno_path}/models'
    n=3000
    
    #download model scores from batches
    modelScores = pd.read_csv(f'{scoresPath}/sklearnModelScoresSectionsTemp.csv')
    
    #feature scores to capture data
    importantFeaturesFile = f'{scoresPath}/importantFeaturesPostShapPostTraining.csv'
    allShapValues = f'{scoresPath}/featureShapValues.all.csv'
    
    ### save empty important features df
    if not os.path.exists(importantFeaturesFile):
        #save empty dataframe on first run
        importantFeaturesShap = pd.DataFrame(columns=['feature','coefs','model'])
        
        print(f'creating important features file : {importantFeaturesFile} ...')
        
        with open(importantFeaturesFile,mode='w',newline='') as f:
            importantFeaturesShap.to_csv(f,index=False)
            f.close()
            
    ### save empty important features df
    if not os.path.exists(allShapValues):
        #save empty dataframe on first run
        allFeaturesShap = pd.DataFrame(columns=['feature','coefs','model'])
        
        print(f'creating important features file : {allShapValues} ...')
        
        with open(allShapValues,mode='w',newline='') as f:
            allFeaturesShap.to_csv(f,index=False)
            f.close()
    
    #get iterations for models with auc > .51
    modelIterations = modelScores.loc[modelScores['auc'] > .51][['data_type','iteration']]
    modelIterations = modelIterations.drop_duplicates().iloc[5:]
    
    for i,data_type in zip(modelIterations['iteration'],modelIterations['data_type']):
        model_str = f'{data_type}_{i}'
        filePath = f'{modelPath}/sklearnGradBoostHistClassifier_{model_str}.pkl'
        try:
            with open(filePath, 'rb') as file:
                clfHGB = pickle.load(open(filePath, 'rb'))
            #clfHGB = pickle.load(model)

            mainDf,sectionPairs = download_data(pheno_path,i,data_type,test_path,epi_file,withdrawal_path)
            y = mainDf['PHENOTYPE'] - 1
            Xmain = mainDf.drop(columns=["PHENOTYPE"])
            print('Xmain array = ',Xmain.shape)
            
            if data_type != 'main':
                Xmain = create_epi_df(Xmain,sectionPairs,combo="sum")
                
            topFeatures,featuresZscores = calculate_plot_shap_values(clfHGB,Xmain,y,i,figPath,data_type)
            #get the featureZscores into a dataframe to merge with dfSnps2
            zscores = pd.DataFrame(data=featuresZscores).reset_index()
            zscores.columns=['feature','shap_zscore']
            
            if topFeatures.empty:
                pass
            else:
                topFeatures = topFeatures.reset_index()
                topFeatures['data_type'] = data_type
                
                with open(importantFeaturesFile,mode='a',newline='') as f:
                    topFeatures.to_csv(f,index=False, header=False)
                    f.close()
            

            featuresZscores = featuresZscores.reset_index()
            featuresZscores['data_type'] = data_type
                
            with open(allShapValues,mode='a',newline='') as f:
                featuresZscores.to_csv(f,index=False, header=False)
                f.close()
            
        except FileNotFoundError:
            print('model file for iteration and data type not found ..',model_str)
            
            

            

        

if __name__ == '__main__':
    
    
    #########   FOR DEBUGGING ######
    pheno = 'myocardialInfarction'
    pheno_path = f'/nfs/scratch/projects/ukbiobank/prsInteractive/results/{pheno}'
    training_path = f'/nfs/scratch/projects/ukbiobank/prsInteractive/results/{pheno}/trainingCombined.raw'
    test_path = f'/nfs/scratch/projects/ukbiobank/prsInteractive/results/{pheno}/testCombined.raw'
    epi_file = f'/nfs/scratch/projects/ukbiobank/prsInteractive/results/{pheno}/epiFiles/trainingCombinedEpi.epi.cc.summary.filtered'
    withdrawal_path=f'/nfs/scratch/projects/ukbiobank/prsInteractive/data/withdrawals.csv'
    
    main(pheno_path,test_path,epi_file,withdrawal_path)