#!/usr/bin/env python3

#clean environmental data (clinical markers such as blood chemistry, counts and lifestyle data such as eating and sleep habits etc)
#clean hla data from raw imputation scores to 0, 1, and 2

import pandas as pd
import os
import sys

def clean_environmental(data_path):
    #download and clean data
    pass

def clean_hla(data_path):
    '''input: hla data in csv format with header = Participant ID, HLA imputation
                                                   float (0-2)
              hla header .txt for column names : ukb_hla_v2.txt

       output: hla variant calls in .csv = [0,1,2] with columns = hla loci
              
    '''
    
    hla = pd.read_csv(f'{data_path}/hla_participant.csv',index_col='Participant ID')
    header = pd.read_csv(f'{data_path}/ukb_hla_v2.txt',sep='\t')
    
    #get it into the correct format
    hla = hla['HLA imputation values'].str.split(",", expand = True)
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
    
    hla5.to_csv(f'{data_path}/participant_hla.csv',index=False)

def main(data_path):
    
    #clean_environmental(data_path)
    clean_hla(data_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creating hla data file...")
    parser.add_argument("--data_folder", help="Path to the input data folder")

    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    data_path = args.data_folder or os.environ.get("DATA_PATH")
    print(f"[PYTHON] Reading from: {data_path}")
    

    if not data_path:
        raise ValueError("You must provide a data path via --data_folder or set the DATA_PATH environment variable.")

            
    main(data_path)