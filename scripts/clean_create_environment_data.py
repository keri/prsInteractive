#!/usr/bin/env python3
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
import os
import argparse
import sys

from helper.download import get_dataset, get_epi_columns, get_columns
from helper.data_wrangling import *

def rank_gene_env_features(geneEnvShapleyFile,threshold):
    '''
    input : df with columns [envGeneticFeature,shap_zscore,env_type,geneticFeature,envFeature,main_E,epistatic]
    
    oupt : df with gene_environment_features ranked with Shapley Z scores > 2
        
    '''
    #download data
    df = pd.read_csv(geneEnvShapleyFile)
    
#   if 'Full' in geneEnvShapleyFile:
#       outputFile = geneEnvShapleyFile.replace('Full','')
#
#   else:
    newFile = geneEnvShapleyFile.split('.')[0]
    outputFile = f'{newFile}Filtered.csv'
        
    #get the largest shap_zscore for duplicated values
    df.sort_values(['envGeneticFeature','shap_zscore'],ascending=False,inplace=True)
    
    #drop duplicates in df
    df.drop_duplicates(subset=['envGeneticFeature'],keep='first',inplace=True)
    
    #get the epistatic interactions
    epiDf = df[df['epistatic'] == 1]
    
    #get the main effects
    mainDf = df[df['main_E'] > 0]
    mainDf.loc[mainDf['envFeature'].isna(),'envFeature'] = mainDf['envGeneticFeature']
    
    importantFeatures = epiDf[(epiDf['shap_zscore'] > threshold) | (epiDf['shap_zscore'] < -threshold) ]
    
    finalFeatures = pd.concat([importantFeatures,mainDf],ignore_index=True)
    
    finalFeatures.to_csv(outputFile,index=False)
    
    return finalFeatures

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
        if len(parts) > 1: #ensure there are not single env features
            first_col = parts[0]
            remaining_cols = ','.join(parts[1:])
            
            combined_columns[epiCombo] = envGeneticDf[first_col] * envGeneticDf[remaining_cols]
        
    # Create DataFrame from dictionary all at once
    combinedDf = pd.DataFrame(combined_columns, index=envGeneticDf.index)
    
    return combinedDf
    

def main(phenoPath,withdrawalPath, trainingPath,testPath,holdoutPath,envDf,hlaDf,importantFeaturesFile,epi_combo,threshold):
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
    
    #check to see if features have been reduced
#   fullFeaturesFile = importantFeaturesFile.split('.')[0]+'Full.csv'
#   if os.path.exists(fullFeaturesFile):
#       print(f'{fullFeaturesFile} exists meaning important GxGxE features have previously been reduced')
#       features = rank_gene_env_features(fullFeaturesFile,threshold=threshold)
#   else:
#       #get ranked GxGxE features
    features = rank_gene_env_features(importantFeaturesFile,threshold)
    
    ################ GET MAIN FEATURES TO STORE SEPARATELY  ############
    
    mainFeatures = features[features['main_E'] > 0]['envFeature'].unique().tolist()
    mainEnvDf = envDf[mainFeatures]
    
    allMainFeatures = features['envFeature'].unique().tolist()
    allEnvDf = envDf[allMainFeatures]
    
    ########### GET GxGxE FEATURES TO PROCESS ###############
    
    #filter the main E features
#   features2 = features[features['main_E'] == 0]
#   
#   #filter features with an epistatic interactions
#   envEpiDf = features2[features2['epistatic'] == 1]
#   epiFeatures = envEpiDf[~envEpiDf['geneticFeature'].isna()]['geneticFeature'].unique().tolist()
#   
#   #remove the epi interactions with HLA region
#   epiGenoFeatures = list(set(epiFeatures) - set(hlaDf.columns.tolist()))
#   
#   envFeatures = envEpiDf['envFeature'].unique().tolist()
#   
#   #GxGxE features comma separated
#   epiEnvFeatures = envEpiDf[~envEpiDf['envFeature'].isna()]['envGeneticFeature'].tolist()
#       
#   
#   ###### PROCESS TRAINING ENVIRONMENT DATA TEST AND HOLDOUT SET ###
#   
#   # DOWNLOAD GENETIC DATA AND MERGE WITH ENV DATA
#   #expand the features into a list and filter redundant features
#   expandedSnps = get_epi_snps(set(epiGenoFeatures))
    
    expandedSnps = []
    
    trainingDf = get_dataset(trainingPath,withdrawalPath,expandedSnps,use_chunking=True)
#   trainingDf = create_epi_df(trainingDf, epiGenoFeatures,combo=epi_combo)
    geneEnvTraining = trainingDf.merge(envDf,left_index=True,right_index=True,how='left')
#   geneEnvTraining = geneEnvTraining.merge(hlaDf,left_index=True,right_index=True,how='left')
    

    testDf = get_dataset(testPath, withdrawalPath, expandedSnps, use_chunking=True)
#   testDf = create_epi_df(testDf, epiGenoFeatures,combo=epi_combo)
    geneEnvTest = testDf.merge(envDf,left_index=True,right_index=True,how='left')
#   geneEnvTest = geneEnvTest.merge(hlaDf,left_index=True,right_index=True,how='left')
    
    
    holdoutDf = get_dataset(holdoutPath, withdrawalPath, expandedSnps,use_chunking=True)
#   holdoutDf = create_epi_df(holdoutDf, epiGenoFeatures,combo=epi_combo)
    geneEnvHoldout = holdoutDf.merge(envDf,left_index=True,right_index=True,how='left')
#   geneEnvHoldout = geneEnvHoldout.merge(hlaDf,left_index=True,right_index=True,how='left')
    
    ######### GET THE MAIN ENV DF #######################
    
    ### Main env with highly protective or risk predictions alone (most predictive in model) or predictive above shap z-score > 2 and in a GxGxE interaction for phenotype  ########
    
    trainingMainDf = mainEnvDf.reindex(trainingDf.index)
    testMainDf = mainEnvDf.reindex(testDf.index)
    holdoutMainDf = mainEnvDf.reindex(holdoutDf.index)

    
    ### env features trained with model ###############
    
    trainingAllEnvDf = allEnvDf.reindex(trainingDf.index)
    testAllEnvDf = allEnvDf.reindex(testDf.index)
    holdoutAllEnvDf = allEnvDf.reindex(holdoutDf.index)

    
#   # MEAN CENTER
#   
#   envTraining,training_mean = mean_center_data(geneEnvTraining[envFeatures])
#   envTest = geneEnvTest[envFeatures] - training_mean
#   envHoldout = geneEnvHoldout[envFeatures] - training_mean
#
#   
#   # REPLACE ENV FEATURES WITH MEAN CENTERED FEATURES
#   geneEnvTraining[envFeatures] = envTraining[envFeatures]
#   geneEnvTest[envFeatures] = envTest[envFeatures]
#   geneEnvHoldout[envFeatures] = envHoldout[envFeatures]
#   
#   # CREATE GENE-ENVIRONMENT FEATURE DATAFRAME
#   combinedTraining = combine_gene_environment(geneEnvTraining,epiEnvFeatures)
#   combinedTest = combine_gene_environment(geneEnvTest,epiEnvFeatures)
#   combinedHoldout = combine_gene_environment(geneEnvHoldout,epiEnvFeatures)
    
    
    # SCALE GENE-ENVIRONEMNT FEATURE DATAFRAME
#   combinedTraining,scaler_model = scale_cardio_training_data(combinedTraining)
#   scaled_test = scaler_model.transform(combinedTest)
#   combinedTest = pd.DataFrame(scaled_test, columns=combinedTest.columns,index=combinedTest.index)
#   scaled_holdout = scaler_model.transform(combinedHoldout)
#   combinedHoldout = pd.DataFrame(scaled_holdout, columns=combinedHoldout.columns,index=combinedHoldout.index)
    
    
    #Main ENV dataframe
    ### MEAN CENTER #########
    mainTraining, mainTraining_mean = mean_center_data(trainingMainDf)
    mainTest = testMainDf[mainFeatures] - mainTraining_mean
    mainHoldout = holdoutMainDf[mainFeatures] - mainTraining_mean

    ####### SCALE DATAFRAME #######
    scaledMainTraining,scaler_model = scale_cardio_training_data(mainTraining)
    scaled_test = scaler_model.transform(mainTest)
    scaledMainTest = pd.DataFrame(scaled_test, columns=scaledMainTraining.columns,index=testMainDf.index)
    scaled_holdout = scaler_model.transform(mainHoldout)
    scaledMainHoldout = pd.DataFrame(scaled_holdout, columns=scaledMainTraining.columns,index=holdoutMainDf.index)
    
    #ALL ENV dataframe MAIN and in EPI interactions
    
    ### MEAN CENTER #########
    allEnvTraining, allEnvTraining_mean = mean_center_data(trainingAllEnvDf)
    allEnvTest = testAllEnvDf[allEnvDf.columns] - allEnvTraining_mean
    allEnvHoldout = holdoutAllEnvDf[allEnvDf.columns] - allEnvTraining_mean

    
    scaledAllEnvTraining,scaler_model = scale_cardio_training_data(allEnvTraining)
    scaled_test = scaler_model.transform(allEnvTest)
    scaledAllEnvTest = pd.DataFrame(scaled_test, columns=scaledAllEnvTraining.columns,index=allEnvTest.index)
    scaled_holdout = scaler_model.transform(allEnvHoldout)
    scaledAllEnvHoldout = pd.DataFrame(scaled_holdout, columns=scaledAllEnvTraining.columns,index=allEnvHoldout.index)
    

#   combinedTraining.reset_index().rename(columns={'index': 'IID'}).to_csv(f'{phenoPath}/geneEnvironmentTraining.csv',index=False)
#   combinedTest.reset_index().rename(columns={'index': 'IID'}).to_csv(f'{phenoPath}/geneEnvironmentTest.csv',index=False)
#   combinedHoldout.reset_index().rename(columns={'index': 'IID'}).to_csv(f'{phenoPath}/geneEnvironmentHoldout.csv',index=False)
    
    scaledMainTraining.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/mainEnvironmentTraining.csv',index=False)
    scaledMainTest.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/mainEnvironmentTest.csv',index=False)
    scaledMainHoldout.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/mainEnvironmentHoldout.csv',index=False)
    
    scaledAllEnvTraining.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/allEnvironmentTraining.csv',index=False)
    scaledAllEnvTest.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/allEnvironmentTest.csv',index=False)
    scaledAllEnvHoldout.reset_index().rename(columns={'index':'IID'}).to_csv(f'{phenoPath}/allEnvironmentHoldout.csv',index=False)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="creating GxGxE Dataframes...")
    parser.add_argument("--training_file", help="data file of training data")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--holdout_file", help="data file of holdout data")
    parser.add_argument("--pheno_data", help="Results Path to write to")
    parser.add_argument("--env_file",help="Environmental data for participants")
    parser.add_argument("--hla_file",help="HLA data for participants")
    parser.add_argument("--gene_env_file",help="Genetic environmental epistatic features")
    parser.add_argument("--withdrawal_path",help="Genetic withdrawal path for IDs")
    parser.add_argument("--threshold", type=float, default=1.99, help="shapley z score threshold to use for filter (default: %(default)s)")
    parser.add_argument("--epi_combo", default="sum", help="env variable describing how epi SNPs are combined (sum, prod)")

    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    training_file = args.training_file or os.environ.get("TRAINING_PATH")
    print(f"training file : {training_file}")
    
    test_file = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_file}")
    
    holdout_file = args.holdout_file or os.environ.get("HOLDOUT_PATH")
    print(f"holdout file : {holdout_file}")
    
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
    env_file = args.env_file or os.environ.get("ENV_FILE")
    print(f"reading from participant environment file : {env_file}")
    
    hla_file = args.hla_file or os.environ.get("HLA_FILE")
    print(f"reading from participant hla file : {hla_file}")
    
    gene_env_file = args.gene_env_file or os.environ.get("GENE_ENV_FILE")
    print(f"reading from gene environment epi features file : {gene_env_file}")
    
    withdrawal_path = args.withdrawal_path or os.environ.get("WITHDRAWAL_PATH")
    print(f"reading withdrawals from file : {withdrawal_path}")
    
    threshold = os.environ.get("THRESHOLD")
    threshold = float(threshold) if threshold else args.threshold
    print(f"analyzing top features based on shap z score : {threshold}")
    
    epi_combo = args.epi_combo or os.environ.get("EPI_COMBO")
    print(f"epi combo to use for combining epi SNPS : {epi_combo}")
#   
    ###########  TEST VARIABLES ##########
    pheno = 'type2Diabetes'
    pheno_data = f"/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi"
    env_file = "/Users/kerimulterer/prsInteractive/results/participant_environment.csv"
    hla_file = "/Users/kerimulterer/prsInteractive/results/participant_hla.csv"
    training_file = f"/Users/kerimulterer/prsInteractive/results/{pheno}/trainingCombined.raw"
    test_file = f"/Users/kerimulterer/prsInteractive/results/{pheno}/testCombined.raw"
    holdout_file = f"/Users/kerimulterer/prsInteractive/results/{pheno}/holdoutCombined.raw"
    gene_env_file=f"{pheno_data}/scores/cardioMetabolicimportantFeaturesPostShap.csv"
    withdrawal_path = '/Users/kerimulterer/prsInteractive/data/withdrawals.csv'
    epi_combo = 'sum'
    threshold=1.99
    
    print(f"[PYTHON] Reading training data from: {training_file}")
    print(f"[PYTHON] Reading test data from: {test_file}")
    print(f"[PYTHON] Reading holdout data from: {holdout_file}")
    print(f"[PYTHON] Writing to phenotype output folder: {pheno_data}")
    print(f"[PYTHON] Reading environmental data from: {env_file}")
    print(f"[PYTHON] Reading HLA data from: {hla_file}")
    print(f"[PYTHON] Reading gene environmental feature data from: {gene_env_file}")
    print(f"Reading threshold for filtering: {threshold}")
    
    
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

        

    
    envDf = pd.read_csv(env_file)
    envDf.rename(columns={'Participant ID':'IID'},inplace=True)
    envDf.set_index(['IID'],inplace=True)

    
    hlaDf = pd.read_csv(hla_file)
    hlaDf.rename(columns={'Participant ID':'IID'},inplace=True)
    hlaDf.set_index(['IID'],inplace=True)
    
    
    main(pheno_data,withdrawal_path,training_file,test_file,holdout_file,envDf,hlaDf,gene_env_file,epi_combo,threshold=threshold)
#   main(args.pheno_data,args.training_file,args.test_file,args.holdout_file,args.env_type,envDf,hlaDf)
        

    
    


    