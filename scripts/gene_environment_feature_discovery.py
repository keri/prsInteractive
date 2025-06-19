#!/usr/bin/env python3
#Keri Multerer Feb 2024 : 

import pandas as pd
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score, matthews_corrcoef, log_loss, jaccard_score, hamming_loss
import matplotlib.pyplot as plt
import warnings
import time
import pickle
import os
import csv
import argparse
from helper.data_wrangling import *
from helper.download import *
from helper.calculate_shap_values import *
#from helper.find_important_features_from_shap import *
warnings.filterwarnings('ignore')

def train_models(X,y,modelPath,pheno,env_type,i):
    '''input : X from training dataset and y
    output : pickled model saved to output file'''
    
    print('training models .....')
    
    ##########################
    # Gradient boosted Hist classifier
    ##########################
    st = time.time()
    parameters_hgb = [{'max_iter':[1000,1500,2000],'learning_rate':[.001,.01,.1,1],'l2_regularization': [0,.1,.5]}]
    clfHGB = HistGradientBoostingClassifier(early_stopping='auto')
    grid_search_hgb = GridSearchCV(estimator=clfHGB,param_grid=parameters_hgb,scoring='roc_auc',cv=3,n_jobs=18)
    #	grid_search_hgb.fit(X.to_numpy(),y.to_numpy())
    grid_search_hgb.fit(X,y)
    pickle.dump(grid_search_hgb, open(f'{modelPath}/sklearnGradBoostHistClassifier_{env_type}_{i}.pkl', 'wb'))
    en = time.time()
    timeHGB = (en-st)/60
    print('time if took to train HGB = ',timeHGB,' minutes')
    
    return(grid_search_hgb)

def score_models(X,y,pheno,env_type,modelFile,i,clfHGB):
    '''load pickled models and score with test set'''
    print('scoring models .....')
    
    st = time.time()
    
    #get the feature names for model
    
    #########################
    # Gradient Boosted Hist 
    #########################

    score = clfHGB.score(X, y)
    yHat = clfHGB.predict(X)
    balanced_score = balanced_accuracy_score(y,yHat)
    auc = roc_auc_score(y, clfHGB.predict_proba(X)[:, 1])
    mcc = matthews_corrcoef(y,yHat)
    logloss = log_loss(y,clfHGB.predict_proba(X)[:, 1])
    jscore = jaccard_score(y,yHat)
    hloss = hamming_loss(y,yHat)
    f1score = f1_score(y,yHat)
    fields=['gradient boosted classifier',score,balanced_score,auc,mcc,logloss,jscore,hloss,f1score,env_type,i]
    
    with open(modelFile,mode='a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        f.close()
        
    en = time.time()
    timeTotal = (en-st)/60
    print('time if took to score models = ',timeTotal,' minutes')
    return(auc)



def main(pheno,env_type,phenoPath,trainingPath,testPath,resultsPath):
    ############################   ENVIRONMENT VARIABLES     ######################

    figPath = f'{phenoPath}/figures'
    scoresPath = f'{phenoPath}/scores'
    modelsPath = f'{phenoPath}/models'
    envPath = f'{resultsPath}/participant_environment.csv'
    hlaPath = f'{resultsPath}/participant_hla.csv'
    modelFile = f'{scoresPath}/{env_type}ModelScores.csv'
    importantFeaturesFile =f'{scoresPath}/{env_type}importantFeaturesPostShap.csv'
    
    if not os.path.exists(f'{modelFile}'):
        # create empty dataframe to capture the scores and snps in each iteration
        models = pd.DataFrame(columns=['model','test score','balanced score','auc','matthews_corrcoef','log_loss','jaccard_score','hamming_loss','f1_score','env_type','iteration'])
        with open(modelFile,mode='w',newline='') as f:
            models.to_csv(f,index=False)
            f.close()
            
    if not os.path.exists(f'{importantFeaturesFile}'):
        # create empty dataframe to capture the scores and GxGxE features for every E feature
        combinedImportFeatures = pd.DataFrame(columns=['envGeneticFeature','shap_zscore','env_type','geneticFeature','envFeature','main_E','epistatic'])
        with open(importantFeaturesFile,mode='w',newline='') as f:
            combinedImportFeatures.to_csv(f,index=False)
            f.close()

            

            
    ########################  DOWNLOAD ENVIRONMENTAL MARKERS   ##################################


    envDf = pd.read_csv(envPath)
            
    envDf.rename(columns={'Participant ID':'IID'},inplace=True)
    envDf.set_index('IID',inplace=True)
    
    #filter data with more than 5% missingness
    # Calculate the percentage of missing values for each column
    missing_percentage = envDf.isna().mean() * 100
    
    # Keep columns with less than 5% missing values
    envDf = envDf.loc[:, missing_percentage < 5]
    
    ########################  DOWNLOAD GENOTYPED DATA FOR FINAL FILTERED FEATURES USED IN FINAL MODEL  ##################################
    
    modelFeatures = pd.read_csv(f'{scoresPath}/importantFeaturesPostShap.csv')
    chunk_size=2000

    #for a large number of features, split the job into chunks of 2K
    n_chunks = modelFeatures.shape[0] // chunk_size

        
    for chunk in range(n_chunks+1):
        start = n_chunks*chunk_size
        stop = start+chunk_size
        modelFeatures2 = get_epi_snps(modelFeatures['feature'].tolist()[start:stop])
        
    
    
        trainingData = get_dataset(trainingPath,modelFeatures2)
        y = trainingData['PHENOTYPE']
        trainingData = create_epi_df(trainingData,modelFeatures['feature'].tolist())
        #
        testData = get_dataset(testPath,modelFeatures2)
        yTest = testData['PHENOTYPE']
        testData = create_epi_df(testData,modelFeatures['feature'].tolist())
        
        ##################### merge HLA data ###########################
        
        #hla data with entire dataset
        hlaData = download_hla_data(hlaPath)
        print('shape of hla data set',hlaData.shape)
        hlaData.set_index('IID',inplace=True)
        #get features for later use
        hlaFeatures = hlaData.columns.tolist()
        
        trainingData = trainingData.merge(hlaData,left_index=True,right_index=True,how='left')
        print('final shape of training dataframe for all features after hla merge = ',trainingData.shape)
        
        testData = testData.merge(hlaData,left_index=True,right_index=True,how='left')
        print('final shape of test dataframe for all features after hla merge = ',testData.shape)
        
        
        ###############################  CREATE TRAINING DATA COMBINATIONS OF 2 CARDIOMETABOLIC ACROSS ALL MODEL FEATURES  FOR EACH MODEL #############################
        
        
        # Create an empty DataFrame to store pair-wise interactions
    
        trainingDataIndex = pd.DataFrame(index=trainingData.index)
        cardioFeaturesTraining = envDf.merge(trainingDataIndex,left_index=True,right_index=True,how='right') #ensures the rows are in the correct order for combining into interaction dataset
        testDataIndex = pd.DataFrame(index=testData.index)
        cardioFeaturesTest = envDf.merge(testDataIndex,left_index=True,right_index=True,how='right')
        
        
        # Loop through each env feature and run models
        # create data set for each model with 1 E feature, interaction matrix (ExG), and G matrix combined
        for i in range(len(envDf.columns)):
            feature1 = envDf.columns[i]
            
            #	for j in range(len(modelFeatures2)):
            # Create pairwise dataset with the pair-wise interaction
            interactionDataset = cardioFeaturesTraining.iloc[:, i].values.reshape(-1,1) * trainingData
            interactionDatasetTest = cardioFeaturesTest.iloc[:, i].values.reshape(-1,1) * testData
            
            #change the column names to have the cardioMetabolic feature in the name
            interactionColNames = [f'{feature1},{x}' for x in trainingData.columns]
            
            #create column names for both training and test data
            interactionDataset.columns = interactionColNames
            interactionDatasetTest.columns = interactionColNames
            interactionDataset.index = trainingData.index
            interactionDatasetTest.index = testData.index
            
            #create column for single cardio metabolic feature without interaction
            interactionDataset[feature1] = cardioFeaturesTraining.iloc[:, i]
            interactionDatasetTest[feature1] = cardioFeaturesTest.iloc[:, i]
            
            #merge two datasets to train model and calculate SHAP
            trainingData2 = interactionDataset.merge(trainingData,left_index=True,right_index=True)
            print(trainingData2.shape)
            testData2 = interactionDatasetTest.merge(testData,left_index=True,right_index=True)
            print(testData2.shape)
    
            clfHGB = train_models(trainingData2,y,modelsPath,pheno,env_type,feature1)
            auc = score_models(testData2, yTest, pheno, env_type, modelFile, feature1, clfHGB)
            #if auc > .51:
                #columns = index of features and shap_valueZscores
            topFeatures,featureZscores = calculate_plot_shap_values(clfHGB,trainingData2,testData2,i,figPath,env_type)
    
            if featureZscores.empty:
                'There was an issue with calculating feature z scores in calculate_shap_values.py script...'
                
            else:
                #'envGeneticFeature','shap_zscore','env_type','geneticFeature','main_E','epistatic','envFeature'
                featureZscores = featureZscores.reset_index()
                featureZscores['env_type'] = env_type
                featureZscores.columns = ['envGeneticFeature','shap_zscore','env_type']
                #remove the feature string when combined with genetic features
                featureZscores['geneticFeature'] = featureZscores['envGeneticFeature'].apply( lambda x : x.replace(f'{feature1},',''))
                #remove the feature string when alone
                featureZscores.loc[featureZscores['envGeneticFeature'] == feature1,'geneticFeature'] = ''
                featureZscores['envFeature'] = ''
                featureZscores['main_E'] = 0
                featureZscores['epistatic'] = 0
                
                featureZscores.sort_values(['shap_zscore'],ascending=False,inplace=True)
                
                # fill in columns for E_feature
                featureZscores.loc[featureZscores['envGeneticFeature'].str.contains(feature1),'envFeature'] = feature1
                
                
                # find index where envGeneticFeature == feature1
                idxFeature =  featureZscores.index[featureZscores['envGeneticFeature'] == feature1][0]
                
                if idxFeature == 0:
                    featureZscores.loc[idxFeature,'main_E'] = 1
    
                else:
                    # if the env feature has a high shap_value
                    featuresGreaterThanE = featureZscores[:idxFeature]
    
                    
                    # check to the to see if the combined GxG is greater than G, GxG alone
                    #get the indices of all features that contain the E feature
                    #this returns a list of indices for which the combined GxE is present
                    idxCombinedFeatures = featuresGreaterThanE.index[featuresGreaterThanE['envGeneticFeature'].str.contains(feature1)].tolist()
                    
                    importantCombinedIndices = []
                    #check to see if the combined shap_value of GxE feature is > G alone
                    for idx in idxCombinedFeatures:
                        #get the genetic feature combined with E
                        combinedFeature = featuresGreaterThanE.iloc[idx]['envGeneticFeature']
                        idxGenetic = featuresGreaterThanE.index[featuresGreaterThanE['envGeneticFeature'] == combinedFeature.replace(f'{feature1},','')]
                        if idxGenetic.empty: #combined feature is ranked higher
                            importantCombinedIndices.append(idx)
        
                        else:
                            if idxGenetic[0] > idx:
                                importantCombinedIndices.append(idx)
                            
        
                            
                    # set epistatic == True for all indices with GxGxE
                    featureZscores.loc[importantCombinedIndices,'epistatic'] = 1
    
                with open(importantFeaturesFile,mode='a',newline='') as f:
                    featureZscores.to_csv(f,index=False,header=False)
                    f.close()
                    #create_plots(shap_df,explainer,shap_values,figPath,feature1,env_type)
            
            
            
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="running models for GxGxE features...")
    parser.add_argument("--pheno_folder", help="Path to the input pheno folder")
    parser.add_argument("--training_file", help="data file of training data")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--env_type", help="data type to analyze")
    parser.add_argument("--pheno", help="Phenotype to analyze")
    parser.add_argument("--results_path", help="data path to results")
    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    pheno_path = args.pheno_folder or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    
    pheno = args.pheno or os.environ.get("PHENO")
    print(f"[PYTHON] Phenotype : {pheno}")
    
    env_type = args.env_type or os.environ.get("ENV_TYPE")
    print(f"data type : {env_type}")
    
    training_path = args.training_file or os.environ.get("TRAINING_PATH")
    print(f"training file : {training_path}")
    
    test_path = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_path}")
    
    results_path = args.results_path or os.environ.get("RESULTS_PATH")
    print(f"results path : {results_path}")
    

    
    
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_folder or set the PHENO_PATH environment variable.")
        
    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
        
    if not env_type:
        raise ValueError("You must provide a env type code via --env_type or set the ENV_TYPE environment variable.")
        
    if not training_path:
        raise ValueError("You must provide a training path via --training_path or set the TRAINING_PATH environment variable.")
        
    if not test_path:
        raise ValueError("You must provide a test path via --test_path or set the TEST_PATH environment variable.")
        
    if not results_path:
        raise ValueError("You must provide a results path via --results_path or set the RESULTS_PATH environment variable.")
        
#   pheno='type2Diabetes_test'
#   pheno_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}'
#   env_type='cardioMetabolic'
#   training_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}/trainingCombined.raw'
#   test_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}/testCombined.raw'
#   results_path='/Users/kerimulterer/prsInteractive/results'
    
    main(pheno,env_type,pheno_path,training_path,test_path,results_path)
    