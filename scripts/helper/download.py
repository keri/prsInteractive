#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import time
import os
from pathlib import Path


def download_hla_data(hlaPath):
    '''download the data change Participant ID to IID and set as index for merging with genotyped data
    
    '''
    hla = pd.read_csv(hlaPath)
    
    hla.rename(columns={'Participant ID':'IID'},inplace=True)
    
    return(hla)

def get_columns(resultsPath):
    """
    Optimized column extraction with better error handling
    """
    try:
        # Read only the header line for column names
        with open(resultsPath, 'r') as f:
            header_line = f.readline().strip()
            
        # Handle different separators
        if '\t' in header_line:
            snpList = header_line.split('\t')
        else:
            snpList = header_line.split()
            
        print(f'Downloaded {len(snpList)} column names')
        return snpList
    
    except Exception as e:
        print(f"Error reading columns: {e}")
        # Fallback to pandas method
        snpList = pd.read_csv(resultsPath, sep=r"\s+", nrows=1).columns.tolist()
        print(f'Downloaded {len(snpList)} columns using fallback method')
        return snpList
    
def get_column_index(columns_to_get, full_columns):
    """
    Get indices of specific columns in the full column list
    """
    # Create a mapping for faster lookup
    col_to_idx = {col: idx for idx, col in enumerate(full_columns)}
    
    # Get indices for requested columns
    indices = []
    missing_cols = []
    
    for col in columns_to_get:
        if col in col_to_idx:
            indices.append(col_to_idx[col])
        else:
            missing_cols.append(col)
            
    if missing_cols:
        print(f"Warning: {len(missing_cols)} columns not found in data")
        if len(missing_cols) <= 10:
            print(f"Missing columns: {missing_cols}")
            
    print(f"Found {len(indices)} of {len(columns_to_get)} requested columns")
    return indices

def get_dataset_optimized(df_pathway, columns_to_get, chunk_processing=False):
    """
    Optimized dataset loading with memory management
    
    Parameters:
    - df_pathway: path to .raw file
    - columns_to_get: list of SNP columns to extract
    - chunk_processing: if True, process in chunks for very large files
    """
    st = time.time()
    
    print(f"Loading dataset from: {df_pathway}")
    print(f"Requesting {len(columns_to_get)} SNP columns")
    
    # Get all column names
    full_columns = get_columns(df_pathway)
    
    # Add essential columns
    essential_cols = ['IID', 'PHENOTYPE']
    columns_to_get = essential_cols + columns_to_get
    
    # Get column indices
    idxColumns = get_column_index(columns_to_get, full_columns)
    
    if len(idxColumns) == 0:
        raise ValueError("No valid columns found to load")
        
    # Adjust column names to match found indices
    actual_columns = [full_columns[i] for i in idxColumns]
    
    # Load withdrawal data
    machinePath = '/'.join(df_pathway.split('/')[:-3])
    print('machine path to find withdrawals under: ',machinePath)
    withdrawal_path = f'{machinePath}/data/withdrawals.csv'
    
    print(f'Loading withdrawals from: {withdrawal_path}')
    
    try:
        withdrawn = pd.read_csv(withdrawal_path, header=None)
        withdrawn_ids = set(withdrawn[0].astype(int))  # Convert to string and use set for faster lookup
        print(f'Found {len(withdrawn_ids)} withdrawn participants')
    except FileNotFoundError:
        print("Warning: No withdrawal file found")
        withdrawn_ids = set()
        
    # Load the main data
    print("Loading main dataset...")
    
    if chunk_processing:
        # For very large files, process in chunks
        df = load_in_chunks(df_pathway, idxColumns, actual_columns, withdrawn_ids)
    else:
        # Standard loading
        try:
            # Use numpy for faster loading
            mainArray = np.genfromtxt(df_pathway, 
                                    delimiter=None,  # Auto-detect delimiter
                                    dtype=float, 
                                    usecols=idxColumns, 
                                    skip_header=1,
                                    filling_values=np.nan,  # Handle missing values
                                    invalid_raise=False)   # Don't crash on bad values
            
            print(f"Loaded array shape: {mainArray.shape}")
            
            # Convert to DataFrame
            df = pd.DataFrame(data=mainArray, columns=actual_columns)
            
        except Exception as e:
            print(f"numpy loading failed: {e}")
            print("Trying pandas fallback...")
            
            # Fallback to pandas
            df = pd.read_csv(df_pathway, 
                           delimiter=r'\s+', 
                           usecols=actual_columns,
                           dtype={'IID': int},  # Keep IID as string
                           low_memory=False)
            
    # Remove withdrawn participants
    if withdrawn_ids:
        print("Removing withdrawn participants...")
        initial_count = len(df)
        df['IID'] = df['IID'].astype(int)  # Ensure string type for comparison
        df = df[~df['IID'].isin(withdrawn_ids)]
        removed_count = initial_count - len(df)
        print(f"Removed {removed_count} withdrawn participants")
        
    # Set index
    df.set_index(['IID'], inplace=True)
    
    # Optimize data types to save memory
    df = optimize_datatypes(df)
    
    en = time.time()
    processing_time = (en - st) / 60
    
    print(f'Dataset loaded successfully:')
    print(f'  Shape: {df.shape}')
    print(f'  Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB')
    print(f'  Processing time: {processing_time:.2f} minutes')

    return df

def load_in_chunks(file_path, col_indices, col_names, withdrawn_ids, chunk_size=10000):
    """
    Load large files in chunks to manage memory
    """
    print(f"Loading in chunks of {chunk_size} rows...")
    
    chunks = []
    
    # Use pandas for chunked reading
    for chunk in pd.read_csv(file_path, 
                            delimiter=r'\s+',
                            usecols=col_names,
                            chunksize=chunk_size,
                            dtype={'IID': int},
                            low_memory=False):
    
        # Remove withdrawn participants from this chunk
        if withdrawn_ids:
            chunk = chunk[~chunk['IID'].isin(withdrawn_ids)]
            
        if len(chunk) > 0:  # Only keep non-empty chunks
            chunks.append(chunk)
    
        print(f"Processed chunk with {len(chunk)} valid rows")
    
    # Combine all chunks
    print("Combining chunks...")
    df = pd.concat(chunks, ignore_index=True)
    print(f"Final memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    return df

def optimize_datatypes(df):
    """
    Optimize data types to reduce memory usage
    """
    print("Optimizing data types...")
    
    # Convert PHENOTYPE to smallest possible integer
    if 'PHENOTYPE' in df.columns:
        df['PHENOTYPE'] = pd.to_numeric(df['PHENOTYPE'], errors='coerce', downcast='integer')
        
    # Convert SNP columns to smallest possible type
    snp_cols = [col for col in df.columns if col not in ['PHENOTYPE']]
    
    for col in snp_cols:
        if df[col].dtype in ['float64', 'int64']:
            # SNP values are typically 0, 1, 2, so int8 is sufficient
            df[col] = pd.to_numeric(df[col], errors='coerce', downcast='integer')
            
    return df

# Enhanced wrapper function
def get_dataset(df_pathway, columns_to_get, use_chunking=False, optimize_memory=True):
    """
    Main function with additional options
    
    Parameters:
    - df_pathway: path to .raw file
    - columns_to_get: list of SNP columns
    - use_chunking: whether to use chunked processing for large files
    - optimize_memory: whether to optimize data types
    """
    
    # Check file size to decide on chunking
    file_size_gb = Path(df_pathway).stat().st_size / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    
    # Auto-enable chunking for very large files
    if file_size_gb > 5 and not use_chunking:
        print("Large file detected - enabling chunked processing")
        use_chunking = True
        
    return get_dataset_optimized(df_pathway, columns_to_get, chunk_processing=use_chunking)

#def get_column_index(snp_list,full_columns):
#   '''input : path/to/merged_allChromosomes.raw
#      columns_to_get : list: [snp1,snp2,....snpN]'''
#   idx_list = []
#   for snp in snp_list:
#       try:
#           idx_list.append(full_columns.index(snp))
#       except ValueError:
#           pass
#   return(idx_list)
#
#def get_dataset(df_pathway,columns_to_get):
#   '''input : epi snps column_list = ['SNP','BEST_SNP','CHR','BEST_CHR']
#   mainfilepath = filepath to raw file
#   output: dataframe space separated .raw, values: 0,1,2 for values, columns: rsID_MA'''
#   st = time.time()
#   
#   columns_to_get = ['IID','PHENOTYPE'] + columns_to_get
##   full_columns = ['FID','IID','PAT','MAT','SEX','PHENOTYPE'] + full_columns
##   idxColumns = get_column_index(columns_to_get,full_columns)
#   
#   #take out the people that have withdrawn from study
#   machinePath = '/'.join(df_pathway.split('/')[:-3])
#   
##   machinePath = '/'.join(df_pathway.split('/')[:-5])
#   print(machinePath)
#   
#   withdrawn = pd.read_csv(f'{machinePath}/data/withdrawals.csv',header=None)
#   print('withdrawals are in path : ',machinePath)
#   
#   with open(df_pathway,'r') as reader:
#       df = pd.read_csv(df_pathway, delimiter='\s+',usecols=columns_to_get)#max_rows=100
#       #df = pd.read_csv(df_pathway, delimiter='\s+',usecols=columns_to_get, nrows=100)
#   en = time.time()
#       
#       
##   df = pd.DataFrame(data=mainArray,columns=columns_to_get)
#   df2 = df[~df['IID'].isin(withdrawn[0])]
#   
#   
#   df2.set_index(['IID'],inplace=True)
#   
#   print(f'time it took to download entire dataset is ',(en-st)/60, ' minutes')
#   return (df2)
#
#
#def get_columns(resultsPath):
#   snpList = pd.read_csv(resultsPath,sep='/s+',nrows=1)
##    print('getting columns...')
##    print('pathway to file for columns = ',resultsPath)
#   
##   snpList = []
##   with open(f'{resultsPath}/merged_allChromosomes.snplist') as f:
##       reader = csv.reader(f,delimiter='\t')
##       for row in reader:
##           if row:  # skip empty rows
##               snpList.append(row[0])
##    df = pd.read_csv(f'{trainingPath}/merged_allChromosomes.snplist',sep='\t',header=None)
##   full_columns = df[0].tolist()
#
#   print('downloaded columns .....')
#   
#   return(snpList)


def get_epi_columns(epi_filepath):
    '''epiFile columns = [CHR, SNP, N_SIG, N_TOT, PROP, BEST_CHISQ, BEST_CHR, BEST_SNP ]
        use the CHR SNP BEST_CHR BEST_SNP'''
    
    epiDf = pd.read_csv(epi_filepath, sep=r'\s+', usecols=['SNP','BEST_SNP'])
    pairList = (epiDf['SNP'] + ',' + epiDf['BEST_SNP']).tolist()
    return (pairList)


if __name__ == "__main__":
    
    phenoPath = '/Users/kerimulterer/prsInteractive/results/type2Diabetes_test'
    trainingPath = '/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/trainingCombined.raw'
    snpList = get_columns(phenoPath)
#   columns_to_get = snpList
#   
    mainArray = get_dataset(trainingPath,snpList, use_chunking=True)