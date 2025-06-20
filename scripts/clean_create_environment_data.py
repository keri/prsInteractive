#!/usr/bin/env python3
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
import os
import argparse
import sys

from helper.download import get_dataset, get_epi_columns, get_columns
from helper.data_wrangling import *

def scale_cardio_training_data(df):
    # Initialize the StandardScaler
    scaler = StandardScaler()
    
    # Fit the scaler to your data (compute the mean and standard deviation)
    scaler.fit(df)
    
    # Transform the data using the fitted scaler
    scaled_data = scaler.transform(df)
    
    # Create a new DataFrame with the scaled data
    scaled_df = pd.DataFrame(scaled_data, columns=df.columns,index=df.index)
    
    return(scaled_df,scaler)

def mean_center_data(df):
    '''mean center data for all features in featureList
        input : df with individuals in training set to be mean centered 
        output : df with name of column preserved but is mean centered for those features in list'''
    mean = df.mean(axis=0)
    dfMeanCentered = df - mean
    
    return(dfMeanCentered,mean)


def impute_data(df):
    '''
    
        input: dataframe env features, genetic features before combining
               index = IID, no phenotype column
        output: imputed Dataframe (training data), imputed_model

       '''
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean').fit(df)
    dfImp = imp_mean.transform(df)
    return(dfImp, imp_mean)

def combine_gene_environment(envGeneticDf,geneticEnvFeatureList):
    # Create all columns at once using a dictionary
    combined_columns = {}
    
    for epiCombo in geneticEnvFeatureList:
        parts = epiCombo.split(',')
        first_col = parts[0]
        remaining_cols = ','.join(parts[1:])
        
        combined_columns[epiCombo] = envGeneticDf[first_col] * envGeneticDf[remaining_cols]
        
    # Create DataFrame from dictionary all at once
    combinedDf = pd.DataFrame(combined_columns, index=envGeneticDf.index)
    
    return combinedDf
    

def main(phenoPath,trainingPath,testPath,holdoutPath,envDf,hlaDf,importantFeaturesFile):
    '''
    input:
        string trainingPath = absolute path to gentoyped training set
        string testPath = absolute path to genotyped test set
        string holdoutPath = absolute path to genotyped holdout set

        envFeaturePath = absolute path to results from _gene_environment_feature_discovery = {env_type}importantFeaturesPostShap.csv

    output: 
        Dataframes : environment-genetic features combined for separate training, test, and holdout
    '''
    
    ############## DOWNLOAD ENVIRONMENTAL DATA #######################
    

    #the G and GxG SNPs post GxGxE feature discovery
    features = pd.read_csv(importantFeaturesFile)
    #filter the main E features
    features2 = features[features['main_E'] == 0]
    
    #filter features with an epistatic interactions
    envEpiDf = features2[features2['epistatic'] == 1]
    epiFeatures = envEpiDf[~envEpiDf['geneticFeature'].isna()]['geneticFeature'].unique().tolist()
    
    #remove the epi interactions with HLA region
    epiGenoFeatures = list(set(epiFeatures) - set(hlaDf.columns.tolist()))
    
    envFeatures = envEpiDf['envFeature'].unique().tolist()
    
    #GxGxE features comma separated
    epiEnvFeatures = envEpiDf[~envEpiDf['envFeature'].isna()]['envGeneticFeature'].tolist()
        
    
    ###### PROCESS TRAINING ENVIRONMENT DATA TEST AND HOLDOUT SET ###
    
    # DOWNLOAD GENETIC DATA AND MERGE WITH ENV DATA
    #expand the features into a list and filter redundant features
    expandedSnps = get_epi_snps(set(epiGenoFeatures))
    
    trainingDf = get_dataset(trainingPath, expandedSnps)
    trainingDf = create_epi_df(trainingDf, epiGenoFeatures)
    geneEnvTraining = trainingDf.merge(envDf,left_index=True,right_index=True,how='left')
    geneEnvTraining = geneEnvTraining.merge(hlaDf,left_index=True,right_index=True,how='left')
    
    
    testDf = get_dataset(testPath, expandedSnps)
    testDf = create_epi_df(testDf, epiGenoFeatures)
    geneEnvTest = testDf.merge(envDf,left_index=True,right_index=True,how='left')
    geneEnvTest = geneEnvTest.merge(hlaDf,left_index=True,right_index=True,how='left')
    
    
    holdoutDf = get_dataset(holdoutPath, expandedSnps)
    holdoutDf = create_epi_df(holdoutDf, epiGenoFeatures)
    geneEnvHoldout = holdoutDf.merge(envDf,left_index=True,right_index=True,how='left')
    geneEnvHoldout = geneEnvHoldout.merge(hlaDf,left_index=True,right_index=True,how='left')
    
    
    # MEAN CENTER
    envTraining,training_mean = mean_center_data(geneEnvTraining[envFeatures])
    envTest = geneEnvTest[envFeatures] - training_mean
    envHoldout = geneEnvHoldout[envFeatures] - training_mean
    
    # REPLACE ENV FEATURES WITH MEAN CENTERED FEATURES
    geneEnvTraining[envFeatures] = envTraining[envFeatures]
    geneEnvTest[envFeatures] = envTest[envFeatures]
    geneEnvHoldout[envFeatures] = envHoldout[envFeatures]
    
    # CREATE GENE-ENVIRONMENT FEATURE DATAFRAME
    combinedTraining = combine_gene_environment(geneEnvTraining,epiEnvFeatures)
    combinedTest = combine_gene_environment(geneEnvTest,epiEnvFeatures)
    combinedHoldout = combine_gene_environment(geneEnvHoldout,epiEnvFeatures)
    
    # SCALE COMBINED DATA
    combinedTraining,scaler_model = scale_cardio_training_data(combinedTraining)
    scaled_test = scaler_model.transform(combinedTest)
    combinedTest = pd.DataFrame(scaled_test, columns=combinedTest.columns,index=combinedTest.index)
    scaled_holdout = scaler_model.transform(combinedHoldout)
    combinedHoldout = pd.DataFrame(scaled_holdout, columns=combinedHoldout.columns,index=combinedHoldout.index)
    
    

    combinedTraining.reset_index().to_csv(f'{phenoPath}/geneEnvironmentTraining.csv',index=False)
    combinedTest.reset_index().to_csv(f'{phenoPath}/geneEnvironmentTest.csv',index=False)
    combinedHoldout.reset_index().to_csv(f'{phenoPath}/geneEnvironmentHoldout.csv',index=False)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="creating GxGxE Dataframes...")
    parser.add_argument("--training_file", help="data file of training data")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--holdout_file", help="data file of holdout data")
    parser.add_argument("--pheno_path", help="Results Path to write to")
    parser.add_argument("--env_file",help="Environmental data for participants")
    parser.add_argument("--hla_file",help="HLA data for participants")
    parser.add_argument("--gene_env_file",help="Genetic environmental epistatic features")
    
    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    training_file = args.training_file or os.environ.get("TRAINING_PATH")
    print(f"training file : {training_file}")
    
    test_file = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_file}")
    
    holdout_file = args.holdout_file or os.environ.get("HOLDOUT_PATH")
    print(f"holdout file : {holdout_file}")
    
    pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    
    env_file = args.env_file or os.environ.get("ENV_FILE")
    print(f"reading from participant environment file : {env_file}")
    
    hla_file = args.hla_file or os.environ.get("HLA_FILE")
    print(f"reading from participant hla file : {hla_file}")
    
    gene_env_file = args.gene_env_file or os.environ.get("GENE_ENV_FILE")
    print(f"reading from gene environment epi features file : {gene_env_file}")
    
    
    
    print(f"[PYTHON] Reading training data from: {training_file}")
    print(f"[PYTHON] Reading test data from: {test_file}")
    print(f"[PYTHON] Reading holdout data from: {holdout_file}")
    print(f"[PYTHON] Writing to phenotype output folder: {pheno_path}")
    print(f"[PYTHON] Reading environmental data from: {env_file}")
    print(f"[PYTHON] Reading HLA data from: {hla_file}")
    print(f"[PYTHON] Reading gene environmental feature data from: {gene_env_file}")

    
    #Check if participant_environment.csv exists
    if not os.path.exists(env_file):
        print(f"ERROR: participant_environment.csv not found at {env_file}")
        print(f"Available files in {env_file}:")
        try:
            for f in os.listdir(env_file):
                print(f"  - {f}")
        except:
            print("  (cannot list directory)")
            sys.exit(1)

        
    ###########  TEST VARIABLES ##########
#   pheno_path = "/Users/kerimulterer/prsInteractive/results/type2Diabetes"
#   env_data_file = "/Users/kerimulterer/prsInteractive/results/participant_environment.csv"
#   hla_data_file = "/Users/kerimulterer/prsInteractive/results/participant_hla.csv"
#   training_file = "/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/trainingCombined.raw"
#   test_file = "/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/testCombined.raw"
#   holdout_file = "/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/holdoutCombined.raw"
#   gene_env_file="/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/scores/cardioMetabolicimportantFeaturesPostShap.csv"
    
    envDf = pd.read_csv(env_file)
    envDf.rename(columns={'Participant ID':'IID'},inplace=True)
    envDf.set_index(['IID'],inplace=True)
    
    hlaDf = pd.read_csv(hla_file)
    hlaDf.rename(columns={'Participant ID':'IID'},inplace=True)
    hlaDf.set_index(['IID'],inplace=True)
    
    
    main(pheno_path,training_file,test_file,holdout_file,envDf,hlaDf,gene_env_file)
#   main(args.pheno_path,args.training_file,args.test_file,args.holdout_file,args.env_type,envDf,hlaDf)
        

    
    


    