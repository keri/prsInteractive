#!/usr/bin/env python3
import pandas as pd
import os
import argparse

from helper.download import get_dataset, get_epi_columns, get_columns
from helper.data_wrangling import *

def mean_center_data(df):
    '''mean center data for all features in featureList
        input : df with individuals in training set to be mean centered 
        output : df with name of column preserved but is mean centered for those features in list'''
    mean = df.mean(axis=0)
    dfMeanCentered = df - mean
    
    return(dfMeanCentered,mean)


def combine_gene_environment(cardioFull,geneticDf,combinedGeneticE):
    pass


def scale_combine_epi_genetic_features(environFull,cardioGeneticFeaturePath,trainingDf):

    envCopy = environFull.copy()
    
    #get the genetic-env features post modelling
    cardioGeneticDf = pd.read_csv(cardioGeneticFeaturePath)
    
    
    envCopy = envCopy[envCopy.index.isin(trainingDf.index)]
    
    envCopyImputed = impute_data(envCopy)
    
    
    
    for epiCombo in epiCardioFeaturesToCombineList:
        XcenteredInteractions[','.join(epiCombo[:])] = Xcentered[epiCombo[0]] * X[epiCombo[1]]
        
    XepiScaled,scalerEpi = scale_cardio_training_data(XcenteredInteractions)
    
    #scale the epiCentered data after combining
    Xcentered,scalerEpiCentered = scale_cardio_training_data(Xcentered)
    
    Xmain = trainingCardioDf[mainCardioFeatures]
    if not Xmain.empty:
        XmainScaled,scalerMain = scale_cardio_training_data(Xmain)
        
    else:
        print("Main cardio DataFrame is empty. Scaling is skipped.")
        scalerMain = None
        
    #save data in dataframe
    #save means for epi features
    df = pd.DataFrame(data=epiMean)
    df.reset_index(inplace=True)
    df.columns = ['feature','mean']
    models = pd.DataFrame({'scalerEpi':scalerEpi,'scalerMain':scalerMain,'scalerEpiCentered':scalerEpiCentered},index=[0])
    df.to_csv(f'{dataPath}/data/cardioTrainingMeans.csv',index=False)
    models.to_csv(f'{dataPath}/data/cardioTrainingModels.csv',index=False)
        
    return(epiMean,scalerEpi,scalerMain,scalerEpiCentered)

def calculate_mean_scaler_for_cardio_training_data(envDf,genoDf,genoEnvEpiFeatures):
    #get full columns to use for index of section columns
    #cardioDf,epiCardioFeaturesToCenter,mainCardioFeatures,epiCardioFeaturesToCombineList,dataPath,get_dataset
    
    print('downloading geno training data for mean and scaler calculation for cardio data ....')
    

        
    genoEpiFeatures = list(set([epi_geno_combo[1] for epi_geno_combo in epiCardioFeaturesToCombineList]))
        
        
    Xcentered,epiMean = mean_center_data(XEpi)
    
    #multiply the mean centered cardio features with geno features for seed epi pairs
    XcenteredInteractions = pd.DataFrame(index=Xcentered.index)
    
#	#get the uncentered data
#	XuncenteredInteractions = pd.DataFrame(cardioDf[epiCardioFeaturesToCenter].index)
    
    for epiCombo in epiCardioFeaturesToCombineList:
        XcenteredInteractions[','.join(epiCombo[:])] = Xcentered[epiCombo[0]] * X[epiCombo[1]]
        
    XepiScaled,scalerEpi = scale_cardio_training_data(XcenteredInteractions)
    
    #scale the epiCentered data after combining
    Xcentered,scalerEpiCentered = scale_cardio_training_data(Xcentered)
    
    Xmain = trainingCardioDf[mainCardioFeatures]
    if not Xmain.empty:
        XmainScaled,scalerMain = scale_cardio_training_data(Xmain)
        
    else:
        print("Main cardio DataFrame is empty. Scaling is skipped.")
        scalerMain = None
        
    #save data in dataframe
    #save means for epi features
    df = pd.DataFrame(data=epiMean)
    df.reset_index(inplace=True)
    df.columns = ['feature','mean']
    models = pd.DataFrame({'scalerEpi':scalerEpi,'scalerMain':scalerMain,'scalerEpiCentered':scalerEpiCentered},index=[0])
    df.to_csv(f'{dataPath}/data/cardioTrainingMeans.csv',index=False)
    models.to_csv(f'{dataPath}/data/cardioTrainingModels.csv',index=False)
        
    return(epiMean,scalerEpi,scalerMain,scalerEpiCentered)



def process_cardio_data(X,cardioDf,cardioDfFull,epiCardioFeaturesToCenter,mainCardioFeatures,epiCardioFeaturesToCombineList,dataPath,full_columns,get_dataset):
    '''transform cardioDf and geno dataset into a combined, mean centered, standardized dataset
    
    input : X dataframe(training or test, genotyped dataset with HLA region
            epiCardioFeaturesToCenter : list of cardio features involved in epistatic interactions that need to be mean centered before combining with geno data
            mainCardioFeatures : list of main cardio features not involved in epi interactions and dont need mean centering
            epiCardioFeaturesToCombineList : zipped list of cardio features and geno features in each epi pair to create new feature for model
            
    return:
            Xtransformed : dataframe combined centered and standardized dataset 
    '''	
    
    mean,scalerEpi,scalerMain,scalerEpiCentered = calculate_mean_scaler_for_cardio_training_data(cardioDfFull,epiCardioFeaturesToCenter,mainCardioFeatures,epiCardioFeaturesToCombineList,dataPath,full_columns,get_dataset)
    
    #mean center the epiCardioFeatures before combining with genotyped data for model training
    cardioMain = cardioDf[mainCardioFeatures]
    cardioEpi = cardioDf[epiCardioFeaturesToCenter]
    
    
    #mean center main	
    Xcentered = cardioEpi-mean
    
    #multiply the mean centered cardio features with geno features for seed epi pairs
    XcenteredInteractions = pd.DataFrame(index=Xcentered.index)
    
#	#get the uncentered data
#	XuncenteredInteractions = pd.DataFrame(cardioDf[epiCardioFeaturesToCenter].index)
    
    for epiCombo in epiCardioFeaturesToCombineList:
        XcenteredInteractions[','.join(epiCombo[:])] = Xcentered[epiCombo[0]] * X[epiCombo[1]]
        
#	for epiCombo in epiCardioFeaturesToCombineList:
#		XuncenteredInteractions[','.join(epiCombo[:])] = cardioDf[epiCombo[0]] * X[epiCombo[1]]
        
    # Transform the data using the fitted scaler
    scaled_data = scalerEpi.transform(XcenteredInteractions)
    
    # Create a new DataFrame with the scaled data
    XepiScaled = pd.DataFrame(scaled_data, columns=XcenteredInteractions.columns,index=XcenteredInteractions.index)
    
    # Transform epi centered data to scaled epi centered data
    scaled_data = scalerEpiCentered.transform(Xcentered)
    XcenteredScaled = pd.DataFrame(scaled_data, columns=Xcentered.columns,index=Xcentered.index)
    
    
    if not cardioMain.empty:
        if scalerMain:
            #transform Xmain
            scaled_data = scalerMain.transform(cardioMain)
            # Create a new DataFrame with the scaled data
            XmainScaled = pd.DataFrame(scaled_data, columns=cardioMain.columns,index=cardioMain.index)
        else:
            print('There was no scaler for main cardio features, in which case there were no cardio main features to scale')
            
    else:
        print("Main cardio DataFrame is empty. Scaling is skipped.")
        XmainScaled = pd.DataFrame()
        
    return(XmainScaled,XepiScaled,XcenteredScaled)


def create_gene_env_training_data(cardioEpiFeatureFile,trainingEnv,trainingPath):
    '''input :
        datapath to files cardioEpiFeatures file, trainingPath to genotyped data
        trainingEnv: dataframe with env features for training data
        
        output:
        combinedGeneticEnvData = Dataframe()
        trained_imputed_model : sklearn.SimpleImputer model
        trained_mean: float
        trained scaler_model = sklearn.Scaler model
        '''

    

def main(resultsPath,cardioEpiFeatureFile,trainingPath,testPath,holdoutPath,env_type):
    
    importantFeaturesFile =f'{resultsPath}/scores/{env_type}importantFeaturesPostShap.csv'
    
    envDf = pd.read_csv(f'{resultsPath}/participant_environment.csv')
    
    #########   GET G AND GXG FEATURES THAT ARE RANKED AND FILTERED TO COMBINE WITH E ##########
    
    trainingID = pd.read_csv(f'{resultsPath}/trainingID.txt',sep=' ')
    testID = pd.read_csv(f'{resultsPath}/testID.txt',sep=' ')
    holdoutID = pd.read_csv(f'{resultsPath}/holdoutID.txt',sep=' ')
    
    trainingEnv = envDf[envDf['IID'].isin(trainingID[0].tolist())]
    trainingEnv.set_index(['IID'],inplace=True)
    
    testEnv = envDf[envDf['IID'].isin(testID[0].tolist())]
    testEnv.set_index(['IID'],inplace=True)
    
    holdoutEnv = envDf[envDf['IID'].isin(holdoutID[0].tolist())]
    holdoutEnv.set_index(['IID'],inplace=True)
    
    #download datasets to combine
    full_columns = get_columns(resultsPath)
    
    #the G and GxG SNPs post GxGxE feature discovery
    features = pd.read_csv(importantFeaturesFile)
    #filter the main E features
    features2 = features[features['main_E'] == 0]
    
    #filter epistatic 
    envEpiDf = features2[features2['epistatic'] == 1]
    epiFeatures = envEpiDf['geneticFeatures'].unique().tolist()
    
    envFeatures = envEpiDf['envFeature'].unique().tolist()
    
    epiEnvFeatures = envEpiDf['envGeneticFeature'].tolist()
    
    #expand the features into a list and filter redundant features
    expandedSnps = get_epi_snps(set(epiFeatures))
    
    trainingDf = get_dataset(trainingDf, expandedSnps)
    
    #create combined GxG and G dataset
    trainingDf = create_epi_df(trainingDf,epiFeatures)
    
    
    
    #process training data and return mean to center, trained imputation model, and trained scaled model
    impute_model, trained_mean, centered_model, geneEnvTrainingData = create_gene_env_training_data(trainingEnv[envFeatures],trainingDf, epiEnvFeatures)
    geneEnvTrainingData.to_csv(f'{resultsPath}/geneEnvironmentTraining.csv')
    
    del geneEnvTrainingData
    
    #use trained models for imputation, centering and scaling validation and holdout sets
    testDf = get_dataset(testPath, expandedSnps)
    imputedDf = impute_model
    
    


    