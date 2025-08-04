#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob

dataPath = '/Users/kerimulterer/ukbiobank'
cardiometabolicDataEpi = pd.read_csv(f'{dataPath}/type2Diabetes/tanigawaSet/featureScores/SHAPFeatureImportance/featureShapValue_cardioMetabolicGxE_epiFeatures.csv',header=None)
cardiometabolicDataEpi.columns = ['pair','shapValue','cardioFeature']
#some of the datasets were run twice so drop duplicates
cardiometabolicDataEpi.drop_duplicates(subset=['pair','cardioFeature'],keep='first',inplace=True)
#create genotyped data column for use in model training when combining filtered cardio-geno pairs
cardiometabolicDataEpi['genoData'] = cardiometabolicDataEpi[['pair','cardioFeature']].apply( lambda row : ','.join(row['pair'].split(',')[len(row['cardioFeature'].split(',')):]),axis=1)
#use only the pairs that don't inlcude the 'Dia BP' feature which is a duplicate of the 'Diastolic blood pressure'
cardioEpiFeatures = cardiometabolicDataEpi.loc[cardiometabolicDataEpi['cardioFeature'] != 'Dia BP']['cardioFeature'].unique()
cardioEpiPairs = cardiometabolicDataEpi.loc[cardiometabolicDataEpi['cardioFeature'] != 'Dia BP']['pair'].unique()

cardiometabolicImportantFeatures = pd.read_csv(f'{dataPath}/type2Diabetes/tanigawaSet/featureScores/SHAPFeatureImportance/importantCardioMetabolicFeaturesPostShap.csv',header=None)
cardiometabolicImportantFeatures.columns = ['pair','cardioFeature','dataType']


finalCardioMainForModeling = cardiometabolicImportantFeatures.loc[(cardiometabolicImportantFeatures['pair'] == cardiometabolicImportantFeatures['cardioFeature'])]['cardioFeature'].unique()
cardioFeaturesToCombine = list(finalCardioMainForModeling) + list(cardioEpiFeatures)
finalFeaturesForModel = list(cardioEpiPairs) + list(finalCardioMainForModeling)

#################################################  END ALL OF US CLEANING  ###############################

#add main effect features to epistatic features to be used in modeling
all = glob.glob(dataPath+'/tanigawaData/metabolomicsBloodChemistry/*_participant.csv') + ['/Users/kerimulterer/ukbiobank/tanigawaData/allOfUsCardioFeatures/AllofUs_participant.cleaned.combined.csv']

cardioFeaturesFinalModel = pd.DataFrame(columns=['cardioFeature','cardioFile'])

combinedDataset = pd.DataFrame()
for f in all:
    columns = pd.read_csv(f,nrows=1).columns.tolist()
    combineFeatures = [feature for feature in cardioFeaturesToCombine if feature in columns]
    if len(combineFeatures) > 0:
        df = pd.read_csv(f,usecols=['Participant ID'] + combineFeatures)
        cardioFeaturesFinalModelTemp = pd.DataFrame({'cardioFeature':combineFeatures,'cardioFile':[f]*len(combineFeatures)})
        cardioFeaturesFinalModel = pd.concat([cardioFeaturesFinalModel,cardioFeaturesFinalModelTemp],ignore_index=True)
        if combinedDataset.empty:
            combinedDataset = df.copy()
        else:
            combinedDataset = combinedDataset.merge(df,on=['Participant ID'])
            
print(combinedDataset.head())
print(combinedDataset.shape)

#get the final cardio feature and cardio feature pairs to be included used in modeling
cardiometabolicImportantFeatures2 = pd.DataFrame({'cardioFeatures':finalFeaturesForModel})
#get the genoData where apropriate
cardiometabolicImportantFeatures2 = cardiometabolicImportantFeatures2.merge(cardiometabolicDataEpi.rename(columns={'pair':'cardioFeatures'})[['cardioFeatures','genoData','cardioFeature']],on=['cardioFeatures'],how='left')
#fill in the cardio feature for data in which the feature wasn't included in the epi dataset
cardiometabolicImportantFeatures2.loc[cardiometabolicImportantFeatures2['cardioFeature'].isna(),'cardioFeature'] = cardiometabolicImportantFeatures2['cardioFeatures']

combinedDataset.to_csv(dataPath+'/type2Diabetes/tanigawaSet/featureScores/cardiometabolicFeatures/combinedImportantFeaturesForModeling_participant.csv',index=False)
cardiometabolicImportantFeatures2.to_csv(dataPath+'/type2Diabetes/tanigawaSet/featureScores/cardiometabolicFeatures/combinedImportantFeaturesForModeling.csv',index=False)



        
    
            


#combine the cardiometabolic features into one dataset to be used in modeling

#build a list of files that have cardiometabolic features manually

    

        