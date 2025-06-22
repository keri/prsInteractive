#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import time


def download_hla_data(hlaPath):
    '''download the data change Participant ID to IID and set as index for merging with genotyped data
    
    '''
    hla = pd.read_csv(hlaPath)
    
    hla.rename(columns={'Participant ID':'IID'},inplace=True)
    
    return(hla)

def get_column_index(snp_list,full_columns):
    '''input : path/to/merged_allChromosomes.raw
       columns_to_get : list: [snp1,snp2,....snpN]'''
    idx_list = []
    for snp in snp_list:
        try:
            idx_list.append(full_columns.index(snp))
        except ValueError:
            pass
    return(idx_list)

def get_dataset(df_pathway,columns_to_get):
    '''input : epi snps column_list = ['SNP','BEST_SNP','CHR','BEST_CHR']
    mainfilepath = filepath to raw file
    output: dataframe space separated .raw, values: 0,1,2 for values, columns: rsID_MA'''
    st = time.time()
    
    columns_to_get = ['IID','PHENOTYPE'] + columns_to_get
#   full_columns = ['FID','IID','PAT','MAT','SEX','PHENOTYPE'] + full_columns
#   idxColumns = get_column_index(columns_to_get,full_columns)
    
    #take out the people that have withdrawn from study
    machinePath = '/'.join(df_pathway.split('/')[:-3])
    
#   machinePath = '/'.join(df_pathway.split('/')[:-5])
    print(machinePath)
    
    withdrawn = pd.read_csv(f'{machinePath}/data/withdrawals.csv',header=None)
    print('withdrawals are in path : ',machinePath)
    
    with open(df_pathway,'r') as reader:
        df = pd.read_csv(df_pathway, delimiter='\s+',usecols=columns_to_get)#max_rows=100
        #df = pd.read_csv(df_pathway, delimiter='\s+',usecols=columns_to_get, nrows=100)
    en = time.time()
        
        
#   df = pd.DataFrame(data=mainArray,columns=columns_to_get)
    df2 = df[~df['IID'].isin(withdrawn[0])]
    
    
    df2.set_index(['IID'],inplace=True)
    
    print(f'time it took to download entire dataset is ',(en-st)/60, ' minutes')
    return (df2)


def get_columns(resultsPath):
    snpList = pd.read_csv(resultsPath,sep='/s+',nrows=1)
#    print('getting columns...')
#    print('pathway to file for columns = ',resultsPath)
    
#   snpList = []
#   with open(f'{resultsPath}/merged_allChromosomes.snplist') as f:
#       reader = csv.reader(f,delimiter='\t')
#       for row in reader:
#           if row:  # skip empty rows
#               snpList.append(row[0])
#    df = pd.read_csv(f'{trainingPath}/merged_allChromosomes.snplist',sep='\t',header=None)
#   full_columns = df[0].tolist()

    print('downloaded columns .....')
    
    return(snpList)


def get_epi_columns(epi_filepath):
    '''epiFile columns = [CHR, SNP, N_SIG, N_TOT, PROP, BEST_CHISQ, BEST_CHR, BEST_SNP ]
        use the CHR SNP BEST_CHR BEST_SNP'''
    
    epiDf = pd.read_csv(epi_filepath, sep='\s+', usecols=['SNP','BEST_SNP'])
    pairList = (epiDf['SNP'] + ',' + epiDf['BEST_SNP']).tolist()
    return (pairList)


if __name__ == "__main__":
    
    phenoPath = '/Users/kerimulterer/prsInteractive/results/type2Diabetes_test'
    trainingPath = '/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/trainingCombined.raw'
    snpList = get_columns(phenoPath)
    columns_to_get = snpList
#   
    mainArray = get_dataset(trainingPath,columns_to_get)