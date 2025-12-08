#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob

from helper.download import *
from helper.calculate_prs import *
from helper.data_wrangling import *
#to draw the scatter plot against prs_main
from helper.graph_grouped_prs_qq_plots import *


def main(pheno):
    #process saige data and calculate prs using holdout set
    
    #columns in saige dataset across 1-22 chromosomes
    #CHR	POS	MarkerID	 Allele1	 Allele2	AC_Allele2	AF_Allele2	MissingRate	BETA	  SE	  Tstat	var	p.value	p.value.NA	Is.SPA	AF_case	AF_ctrl	N_case	N_ctrl	N_case_hom	N_case_het	N_ctrl_hom	N_ctrl_het
    
    #columns we will use: MarkerID,p.value
    
    ################  PROCESS AND COMBINE SUMMARY STATISTICS FILES INTO ONE AND FILTER SNPs BASED ON 10-8 PVALUE  #########

    machinePath = '/Users/kerimulterer'
    dfPathway = f'{machinePath}/ukbiobank/{pheno}/tanigawaSet'
    testPathway = f'{dfPathway}/data/testCombined_final.raw'
    figurePath = f'{dfPathway}/figures/validation/otherMethods'
    
    resultsPath = f'{dfPathway}/prs/otherMethods'
    
    saige = pd.read_csv(f'{resultsPath}/combinedPRSValidationOtherMethods.csv')
    
    meanDiff = saige.groupby(['PHENOTYPE']).mean()['scaled_prs_saige'].diff().loc[2]
    
    
    files =  glob.glob(f'{machinePath}/ukbiobank/{pheno}/results/saige/*.txt')
    
    #dont include the index file output
    files = [f for f in files if 'index' not in f]
    
    #dont include the variance file output
    files = [f for f in files if 'varianceRatio' not in f]
    

    #create empty dataframe to hold all SNPs that pass filtering
    summary_stats = pd.DataFrame()
    
    for f in files:
        df = pd.read_csv(f,sep="\t",usecols=['CHR','MarkerID','BETA', 'SE','p.value'])
        print(df.shape)
        df2 = df[df['p.value'] < .00000005]
        print(df2.shape)
        summary_stats = pd.concat([summary_stats,df2],ignore_index=True)
        
    print(summary_stats.shape)
    
    summary_stats.to_csv(f'{resultsPath}/saige.filtered.csv',index=False)
    ####################  CALCULATE PRS FOR TEST DATA ################
    
    #get full columns to use for index of section columns
    full_columns = get_columns(dfPathway)
    
    
    #CREATE A LIST OF COLUMNS TO DOWNLOAD
    separatedSNPs = summary_stats['MarkerID'].tolist()
    
    mainEpiDf = get_dataset(testPathway,separatedSNPs,full_columns)
    
    #create feature score dictionary to be used for generating prs score
    #takes 2 columns [feature,lasso_coefs (unfortunate name used as all beta columns need to be changed to lasso_coefs to calculate PRS using this function)]
    summary_stats.rename(columns={'MarkerID':'feature','BETA':'lasso_coefs'},inplace=True)
    
    #create feature score dictionary
    featureScoreDict = get_feature_coef_dictionaries(summary_stats)
    
    calculate_create_prs_plots(mainEpiDf,featureScoreDict,'saige.analysis',figurePath,resultsPath,'mixed',separatedSNPs)
    combine_with_other(pheno,resultsPath)
    
    
def combine_with_other(pheno,resultsPath):
    '''calculate the odds ratio and prevalence for saige and combine with previous validation methods
        input: pheno = string 
               resultsPath = file path to other methods results folder
    '''
    try:
        #function saves the prs to resultsPath with columns [PHENOTYPE,IID,prs,scaled_prs + all features]
        prsSaige = pd.read_csv(f'{resultsPath}/saige.analysis.mixed.prs.csv',usecols=['IID','prs','scaled_prs','PHENOTYPE'])
    except FileNotFoundError:
        print("Need to run the prs analysis for saige. Will be few minutes ....")
        #run the main function
        main(pheno,resultsPath)
    
    #bin prs 
    prsSaige['bin_saige'] = pd.qcut(prsSaige['scaled_prs'], 10, labels=list(range(1,11)),duplicates='drop')
    prsSaige.rename(columns={'prs':'prs_saige','scaled_prs':'scaled_prs_saige'},inplace=True)
    #calculate odds ratio top 20 v rest
    odds_ratio = calculate_odds_ratio_for_prs(prsSaige, 'scaled_prs_saige', prscr=False)
    odds_ratio['model'] = 'saige'
    #calculate prevalence
    prevalence = calculate_prevalance(prsSaige, 'bin_saige','PHENOTYPE')
    prevalence['model'] = 'saige'
    prevalence['method'] = 'saige'
    
    try:
        prsCombined = pd.read_csv(f'{resultsPath}/combinedPRSValidationOtherMethods.csv')
        prevalenceDf = pd.read_csv(f'{resultsPath}/prevalenceValidationOtherMethods.csv')
        odds_ratio_df = pd.read_csv(f'{resultsPath}/oddsRatiosValidationOtherMethods.csv')
        
        #check for previous saige columns in dataframe
        saige_columns = [col for col in prsCombined.columns if 'saige' in col]
        if len(saige_columns) > 0:
            prsCombined.drop(columns=saige_columns,inplace=True)
            prevalenceDf = prevalenceDf[~prevalenceDf['method'].str.contains('saige')]
            odds_ratio_df = odds_ratio_df[~odds_ratio_df['model'].str.contains('saige')]
            
        

                        
        #concatenate new values to odds ratio and prevalence
        prsCombined = prsCombined.merge(prsSaige,on=['IID','PHENOTYPE'])
        odds_ratio_df = pd.concat([odds_ratio,odds_ratio_df],ignore_index=True)
        prevalenceDf = pd.concat([prevalenceDf,prevalence],ignore_index=True)
            
        prsCombined.to_csv(f'{resultsPath}/combinedPRSValidationOtherMethods.csv',index=False)
        prevalenceDf.to_csv(f'{resultsPath}/prevalenceValidationOtherMethods.csv',index=False)
        odds_ratio_df.to_csv(f'{resultsPath}/oddsRatiosValidationOtherMethods.csv',index=False)
            
    except FileNotFoundError:
        print("Need to run the pgs prs or file is not in correct location or misspelled")
            

if __name__ == '__main__':
    
#   pheno = sys.argv[1]
    pheno = 'type2Diabetes'
    
    main(pheno)
    
    