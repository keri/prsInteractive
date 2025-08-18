import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import time
import warnings
warnings.simplefilter(action='ignore')
import sys
import os
import argparse
# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)
from helper.download import *
from helper.calculate_prs import *
from helper.draw_plots import *
from helper.data_wrangling import *



def create_prs_direction(df,featureScores,image_str,figurePath,prsPath):
    dfCopy = df.copy()
    featureScoreDict = get_feature_coef_dictionaries(featureScores)
    direction='mixed'
#   for direction in ['risk','protect','mixed']:
#
#       if direction == 'protect':
#           featureScoresTemp = featureScores[featureScores['coefs'] < 0]
#           features = featureScoresTemp['feature'].tolist()
##           featureScoreDict = get_feature_coef_dictionaries(featureScoresTemp)
#           print('# of features in protect PRS :',len(features))
#           meanDiff = calculate_create_prs_plots(dfCopy,featureScoreDict,image_str,figurePath,prsPath,direction,features)
#       elif direction == 'risk':
#           featureScoresTemp = featureScores[featureScores['coefs'] > 0]
#           features = featureScoresTemp['feature'].tolist()
#           print('# of features in risk PRS :',len(features))
#           
##           featureScoreDict = get_feature_coef_dictionaries(featureScoresTemp)
#           meanDiff = calculate_create_prs_plots(dfCopy,featureScoreDict,image_str,figurePath,prsPath,direction,features)
#       else:
#           featureScoresTemp = featureScores[featureScores['coefs'] != 0]
#           #mixed for both main and epi
#           features = featureScoresTemp['feature'].tolist()
#           print('# of features in mixed PRS :',len(features))
#           #sort featureScores
##           featureScoreDict = get_feature_coef_dictionaries(featureScoresTemp)
#           meanDiff = calculate_create_prs_plots(dfCopy,featureScoreDict,image_str,figurePath,prsPath,direction,features)
    featureScoresTemp = featureScores[featureScores['coefs'] != 0]
    #mixed for both main and epi
    features = featureScoresTemp['feature'].tolist()
    print('# of features in mixed PRS :',len(features))
    #sort featureScores
#           featureScoreDict = get_feature_coef_dictionaries(featureScoresTemp)
    meanDiff = calculate_create_prs_plots(dfCopy,featureScoreDict,image_str,figurePath,prsPath,direction,features)
    featureScoreDict.clear()
    


def create_saturation_plots(df,featureScores,data_type,figurePath,prsPath):
    dfCopy = df.copy()
    featureScores.sort_values(['coefs'],ascending=False,inplace=True)
#   figurePathSaturation = f'{figurePath}/'
#   prsPathSaturation = f'{prsPath}/'
    featureScoreDict = get_feature_coef_dictionaries(featureScores)
    if 'holdout' in data_type:
        scale_prs=False
    else:
        scale_prs=True
    #create saturation plot for increments of 20 if main and 10 if epi
#   for data_type in ['main','epi']:
    if 'main' in data_type:
        n = 20
    elif 'all' in data_type:
        n=20
    else:
        n = 10
#   featureScores2 = featureScores[featureScores['model'] == data_type]
    #protect dataset
    featureScoresRisk = featureScores[featureScores['coefs'] > 0]
    topN = featureScoresRisk.shape[0]
    #risk dataset
    featureScoresProtect = featureScores[featureScores['coefs'] < 0]
    if featureScoresProtect.shape[0] > topN:
        topN = featureScoresProtect.shape[0]
    meanDiffs = []
    directionList = []
    nFeatures = []
    
    #create increment list for loop
    incList = list(range(n,topN,n))
    
    for i in incList:
        #get the scores for protect, risk, and mixed
        #risk
        direction='risk'
        featuresRisk = featureScoresRisk.head(i)['feature'].tolist()
        prsDf, meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresRisk,scale_prs=scale_prs)
        meanDiffs.append(meanDiff)
        directionList.append(direction)
        nFeatures.append(i)
        
        #protect
        direction='protect'
        featuresProtect = featureScoresProtect.tail(i)['feature'].tolist()
        prsDf, meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresProtect,scale_prs=scale_prs)
        meanDiffs.append(meanDiff)
        directionList.append(direction)
        nFeatures.append(i)
        
        #mixed
        direction='mixed'
        featuresBoth = featuresRisk + featuresProtect
        prsDf,meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresBoth,scale_prs=scale_prs)
        meanDiffs.append(meanDiff)
        directionList.append(direction)
        nFeatures.append(i)        
    
        
    #get the last iteration

    #risk
    direction='risk'
    featuresRisk = featureScoresRisk['feature'].tolist()
    prsDf,meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresRisk,scale_prs=scale_prs)
    meanDiffs.append(meanDiff)
    directionList.append(direction)
    nFeatures.append(topN)
    #protect
    direction='protect'
    featuresProtect = featureScoresProtect['feature'].tolist()
    prsDf,meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresProtect,scale_prs=scale_prs)
    meanDiffs.append(meanDiff)
    directionList.append(direction)
    nFeatures.append(topN)
    #mixed
    direction='mixed'
    featuresBoth = featureScores['feature'].tolist()
    prsDf,meanDiff = calculate_prs(dfCopy,featureScoreDict,featuresBoth,scale_prs=scale_prs)
    meanDiffs.append(meanDiff)
    directionList.append(direction)
    nFeatures.append(topN)
    
    #make dataframe to use for plotting
    
    diffDf = pd.DataFrame({'number_features':nFeatures,'risk_direction':directionList,'mean_diff':meanDiffs})
    diffDf.set_index('number_features', inplace=True)
    diffDf.to_csv(f'{prsPath}/{data_type}.saturationPrs.csv')
    #set title 
    title = f'Saturation Plot for {data_type} in batches of {n} features'
    diffDf.groupby('risk_direction')['mean_diff'].plot(figsize=(10,10),legend=True,title=title,xticks=range(0,topN,n),rot=45)
    plt.savefig(f'{figurePath}/{data_type}.saturationPlot.png')
    plt.clf()

        
    featureScoreDict.clear()
    


#(pheno,pheno_path,test_path,test_env_file,holdout_path,holdout_env_file,scores_path,covar_file,hla_file,feature_file,topN=10000)
def main(pheno,withdrawalPath,phenoPath,testPathway,envFileTest,holdoutPathway,envFileHoldout,covarFile,hlaFile,featureFile,topN=10000):
#   pheno = 'type2Diabetes'
#   topN=1000
        ##########################################
        #          CREATE VARIABLES              #
        ##########################################

    
    prsPath = f'{phenoPath}/scores'
    figurePath = f'{phenoPath}/figures'

    ##############################################################################
    #                       DOWNLOAD COVARIATE DATA                              #
    ##############################################################################
    
    #download covariate data
    covDf = download_covar_data(covarFile)
    covDf.set_index('IID',inplace=True)
    
    #NOT TO INCLUDE IN THE PRS/iPRS CALCULATION
    covarFeatures = covDf.columns.tolist()
    
    ############################  GET FEATURES FROM FINAL MODEL THAT HAVE BEEN PRUNED FOR LD  #################
    
    #for the dataset
    filteredFeatures = pd.read_csv(featureFile)    
    
    filteredFeatures = filteredFeatures[~filteredFeatures['feature'].str.contains('Intercept')]
    

    #filter covariate features out of dataset
    covarCoefs = filteredFeatures[filteredFeatures['model'] == 'covariate']
    covarFeatures = covarCoefs['feature'].tolist()
    
    #get all other models that are not covariate
    filteredFeatures = filteredFeatures[filteredFeatures['model'] != 'covariate']

    
    #filters the covar features out of each model type
    filteredFeatures = filteredFeatures[~filteredFeatures['feature'].isin(covarFeatures)]

    #get the main e features
    cardioMainData = filteredFeatures[filteredFeatures['model'] == 'cardio_main']
    cardioMainPlusAll = filteredFeatures[filteredFeatures['model'] == 'all+cardio_main']
    
    #remove main clinical features from data
    filteredFeatures = filteredFeatures[filteredFeatures['model'] != 'cardio_main']
    filteredFeatures = filteredFeatures[filteredFeatures['model'] != 'all+cardio_main']

    #get full columns to use for index of section columns
    full_columns = get_columns(testPathway)
    

    ################################################################################
    #                      DOWNLOAD TEST/HOLDOUT DATASET                           #
    ################################################################################
        
    #get the cardio features which may be different from some of the epi and main features
    allColumns = list(set(filteredFeatures['feature'].tolist()))
    print('total geno features to download = ',len(allColumns))
    
    #separated features to combine in final dataset for PRS calculations. 
    #will include E features as well
    separatedSNPs = get_epi_snps(allColumns)
    print('total number of SNPs to download = ',len(separatedSNPs))
    
    #get the G and GxG only features for download of genotyped data
    geneticFeatures = list(set(filteredFeatures[filteredFeatures['model'].isin(['main','epi'])]['feature'].tolist()))
    featuresToDownload = get_epi_snps(geneticFeatures)
    
    #######################################################
    #                 DOWNLOAD HLA DATA                   #
    #######################################################
    
    #hla data with entire dataset
    hlaData = download_hla_data(hlaFile)
    print('shape of hla data set',hlaData.shape)
    hlaData.set_index('IID',inplace=True)
    #get features for later use
    hlaFeatures = hlaData.columns.tolist()
    
    #cardioEnv data
    cardioEnvTest = pd.read_csv(envFileTest)
    cardioEnvTest.set_index('IID',inplace=True)
    
    cardioEnvHoldout = pd.read_csv(envFileHoldout)
    cardioEnvHoldout.set_index('IID',inplace=True)
    
    

    for holdout_str in ['','.holdout']:
#   for holdout_str in ['']:
            
        if holdout_str == '':
            data_pathway = testPathway
            cardioEnvData = cardioEnvTest
        else:
            data_pathway = holdoutPathway
            cardioEnvData = cardioEnvHoldout
            
        ####################  DOWNLOAD THE DATASET WITH ALL GENO FEATURES ###################
        mainEpiDf = get_dataset(data_pathway,withdrawalPath,featuresToDownload,use_chunking=True)
        
        #get the IIDs for dataset to be used in cardio data and covar PRS
        #get the participants in the mainEpiDf dataset
        indexMatch = pd.DataFrame(index=mainEpiDf.index)
        
        #save to merge again after epi creation
        phenotype = mainEpiDf[['PHENOTYPE']]
        
        #merge hla features with genofeatures 
        mainEpiDf = mainEpiDf.merge(hlaData,left_index=True,right_index=True,how='left')
        print('final shape of training dataframe for all features after hla merge = ',mainEpiDf.shape)
        
        #create combined feature dataframe  
        mainEpiDf = create_epi_df(mainEpiDf,geneticFeatures)
        
        #merge back to the PHENOTYPE column
        mainEpiDf = mainEpiDf.merge(phenotype,left_index=True,right_index=True)
        print('shape of main dataset = ',mainEpiDf.shape)
        
        mainEpiDf = mainEpiDf.merge(cardioEnvData,left_index=True,right_index=True)
        print('shape of main dataset after environment data concatenated = ',mainEpiDf.shape)
        
        #use this to calculate PRS for covar features
        covarDf = phenotype.merge(covDf,left_index=True,right_index=True,how='left')
        print('covariate only df shape = ',covarDf.shape)
        
        
        
        ####################################################
        #           GET COMBINED MAIN AND EPI DATASET      #
        #####################################################
        
        for model in filteredFeatures['model'].unique():
#       for model in ['main','epi','epi+main','cardio']:
#       for model in prs_models:
            #get the feature coefs for each model in turn. The covariate features have been filtered out
            prsFeatures = filteredFeatures[filteredFeatures['model'] == model]
            modelN = prsFeatures.shape[0]
            image_str = f'{model}.{modelN}{holdout_str}' 
            print(f'# of features in {model} model = ',modelN)
            print('image string = ',image_str)
            print('')
            
#           if model == 'cardio_main':
#               create_prs_direction(cardioMainData,prsFeatures,image_str,figurePath,prsPath)
#
#           else:
            if model != "covariate":
                create_prs_direction(mainEpiDf,prsFeatures,image_str,figurePath,prsPath)
                create_saturation_plots(mainEpiDf, prsFeatures, image_str, figurePath, prsPath)
                    
            else:
                create_prs_direction(covarDf,covarCoefs,image_str,figurePath,prsPath)
                
                
            #if model == 'all' then separate epi,main,epi+main,and cardio
            if model == 'all':
                for sub_model in ['epi','main','epi+main','cardio']:
                    subFeatures = filteredFeatures[filteredFeatures['model'] == sub_model]['feature'].tolist()
                    sub_data = prsFeatures[prsFeatures['feature'].isin(subFeatures)]
                    modelN = sub_data.shape[0]
                    image_str = f'{model}.{modelN}{holdout_str}' 
                    create_prs_direction(mainEpiDf,sub_data,f'{image_str}.{sub_model}.FromAll',figurePath,prsPath)
                
        
        #create PRS for covariate only model
#       modelN = covarCoefs.shape[0]
#       image_str = f'covar.{modelN}{holdout_str}'
#       print('')
#       print('image string = ',image_str)
#       print(f'number of features in covar prs ',modelN)
#       create_prs_direction(covarDf,covarCoefs,image_str,figurePath,prsPath)
    



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="calculating PRS for association modelling weights....")
    parser.add_argument("--pheno_folder", help="Path to the input pheno folder")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--holdout_file", help="data file of test data")
    parser.add_argument("--covar_file", help="data file of covar data")
    parser.add_argument("--hla_file", help="data file of hla data")
    parser.add_argument("--test_env_gen_file", help="test environmental data to use")
    parser.add_argument("--holdout_env_gen_file", help="holdout environmental data to use")
    parser.add_argument("--pheno", help="Phenotype to analyze")
    parser.add_argument("--feature_scores_file", help="data path to feature scores from association modelling")
    parser.add_argument("--withdrawal_path",help="Genetic withdrawal path for IDs")
    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    pheno_path = args.pheno_folder or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    
    test_path = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_path}")
    
    holdout_path = args.holdout_file or os.environ.get("HOLDOUT_PATH")
    print(f"holdout file : {holdout_path}")
    
    covar_file = args.covar_file or os.environ.get("COVAR_FILE")
    print(f"covariate file : {covar_file}")
    
    hla_file = args.hla_file or os.environ.get("HLA_FILE")
    print(f"HLA file : {hla_file}")
    
    test_env_file = args.test_env_gen_file or os.environ.get("GENE_ENV_TEST")
    print(f"environmental test data : {test_env_file}")
    
    holdout_env_file = args.holdout_env_gen_file or os.environ.get("GENE_ENV_HOLDOUT")
    print(f"environmental holdout data : {holdout_env_file}")
    
    pheno = args.pheno or os.environ.get("PHENO")
    print(f"[PYTHON] Phenotype : {pheno}")
    
    feature_file = args.feature_scores_file or os.environ.get("FEATURE_SCORES_FILE")
    print(f"feature scores file : {feature_file}")
    
    withdrawal_path = args.withdrawal_path or os.environ.get("WITHDRAWAL_PATH")
    print(f"reading withdrawals from file : {withdrawal_path}")
    

    
#   pheno='type2Diabetes'
#   pheno_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}'
#   test_env_file=f'/Users/kerimulterer/prsInteractive/results/{pheno}/geneEnvironmentTest.csv'
#   holdout_env_file=f'/Users/kerimulterer/prsInteractive/results/{pheno}/geneEnvironmentHoldout.csv'
#   test_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}/testCombined.raw'
#   holdout_path=f'/Users/kerimulterer/prsInteractive/results/{pheno}/holdoutCombined.raw'
#   results_path='/Users/kerimulterer/prsInteractive/results'
#   covar_file='/Users/kerimulterer/prsInteractive/results/covar.csv'
#   hla_file='/Users/kerimulterer/prsInteractive/results/participant_hla.csv'
#   feature_file=f'/Users/kerimulterer/prsInteractive/results/{pheno}/scores/importantFeaturesForAssociationAnalysis.csv'
#   withdrawal_path = f'/Users/kerimulterer/prsInteractive/data/withdrawals.csv'
    
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_folder or set the PHENO_PATH environment variable.")
        
    if not test_path:
        raise ValueError("You must provide a test path via --test_path or set the TEST_PATH environment variable.")

    if not holdout_path:
        raise ValueError("You must provide a holdout path via --test_path or set the HOLDOUT_PATH environment variable.")
        
    if not covar_file:
        raise ValueError("You must provide a covar file path via --covar_file or set the COVAR_FILE environment variable.")

    if not hla_file:
        raise ValueError("You must provide a hla file path via --hla_file or set the HLA_FILE environment variable.")

    if not test_env_file:
        raise ValueError("You must provide a test env file via --test_env_file or set the GENE_ENV_TEST environment variable.")
        
    if not holdout_env_file:
        raise ValueError("You must provide a holdout env data file via --holdout_env_file or set the GENE_ENV_HOLDOUT environment variable.")

    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
        
    if not feature_file:
        raise ValueError("You must provide a feature scores file from association modelling --feature_scores_file or set the FEATURE_SCORES_FILE environment variable.")
    
    if not withdrawal_path:
        raise ValueError("You must provide a path to withdrawals --withdrawal_path or set the WITHDRAWAL_PATH environment variable.")
            
    main(pheno,withdrawal_path,pheno_path,test_path,test_env_file,holdout_path,holdout_env_file,covar_file,hla_file,feature_file,topN=10000)


    










    