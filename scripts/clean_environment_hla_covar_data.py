#!/usr/bin/env python3

#clean environmental data (clinical markers such as blood chemistry, counts and lifestyle data such as eating and sleep habits etc)
#clean hla data from raw imputation scores to 0, 1, and 2

import pandas as pd
import os
import argparse


def rename_covar_columns(covarFile):
    '''rename columns to simple headings'''
    columns = covarFile.columns
    columns1 = [('PC'+ col.split(' ')[-1]) if 'principal' in col else col for col in columns ]
    covarFile.columns = columns1
    covarFile.rename(columns={'Age at recruitment':'age'},inplace=True)
    return(covarFile)

def clean_environmental(data_path):
    #download and clean data
    df = pd.read_csv(f'{data_path}/participant_environment.csv')
    return(df)
    
def create_covar_data(data_path):
    df = pd.read_csv(f'{data_path}/participant.csv')
    df.rename(columns={'Participant ID':'IID'},inplace=True)
    df['SEX'] = df['Sex'].apply( lambda x : 0 if x == 'Female' else 1)
    return(df)
    
    
def clean_hla(data_path):
    '''input: hla data in csv format with header = Participant ID, HLA imputation
                                                   float (0-2)
              hla header .txt for column names : ukb_hla_v2.txt

       output: hla variant calls in .csv = [0,1,2] with columns = hla loci
              
    '''
    
    hla = pd.read_csv(f'{data_path}/hla_participant.csv',index_col='Participant ID')
    header = pd.read_csv(f'{data_path}/ukb_hla_v2.txt',sep='\t')
    
    try:
      #get it into the correct format. Will fail in test run
      hla = hla['HLA imputation values'].str.split(",", expand = True)
    
    except KeyError:
      pass
      
      
    hla = hla.astype(float)
    
    
    #set copy # to 1 with a threshold of .7
    hla1 = hla.mask((hla >= .7) & (hla < 1), 1)
    #set copy # to 0 if threshold < .7
    hla2 = hla1.mask((hla1 < .7), 0)
    #set copy # to 2 with threshold > 1.7 and < 2
    hla3 = hla2.mask((hla2 >= 1.7) & (hla < 2), 2)
    #set copy number to 1 with threshold < 1.7
    hla4 = hla3.mask((hla3 < 1.7) & (hla > 1), 1)
    
    hla4.columns = header.columns
    
    #filter alleles with prevelance less than .1%
    hlaPrevalence = hla4[hla4 > 0].count() / hla4.count()
    hlaPrevalenceFinal = hlaPrevalence[hlaPrevalence > .001]
    hlaToKeep = hlaPrevalenceFinal.index.tolist()
    
    hla5 = hla4[hlaToKeep]
    hla5.reset_index(inplace=True)
    return(hla5)
    

def main(data_path,results_path):
    
    #clean_environmental(data_path)
    hla_data = clean_hla(data_path)
    hla_data.to_csv(f'{results_path}/participant_hla.csv',index=False)
    
    covar_data = create_covar_data(data_path)
    covar_data.to_csv(f'{results_path}/covar.txt', sep=' ', index=False)
    
    env_data = clean_environmental(data_path)
    env_data.to_csv(f'{results_path}/participant_environment.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creating hla and covar data file...")
    parser.add_argument("--data_folder", help="Path to the input data folder")
    
    parser.add_argument("--results_folder", help="Path to the results data folder")

    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    data_path = args.data_folder or os.environ.get("DATA_PATH")

    if not data_path:
        raise ValueError("You must provide a data path via --data_folder or set the DATA_PATH environment variable.")
        
    print(f"[PYTHON] Reading from: {data_path}")
    
    
    results_path = args.data_folder or os.environ.get("RESULTS_PATH")
    
    if not results_path:
        raise ValueError("You must provide a results path via --results_path or set the RESULTS_PATH environment variable.")
    
    print(f"[PYTHON] Writing to: {results_path}")

            
    main(data_path,results_path)