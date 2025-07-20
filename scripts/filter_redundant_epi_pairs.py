#!/usr/bin/env python3
import pandas as pd
import argparse
import os
from helper.download import get_epi_columns

def find_matching_epi_features_but_backwards(df):
    dfCopy = df.copy()
    
    dfCopy['duplicated'] = 0
    for pair in dfCopy['feature'].tolist():
        pairList = pair.split(',')
        new_pair = pairList[1]+','+pairList[0]
        # locate pairs that are present in reverse
        if dfCopy[dfCopy['feature'] == new_pair].shape[0] > 0:
            dfCopy.loc[dfCopy['feature'] == pair,'duplicated'] = 1

            
    filteredDf = dfCopy[dfCopy['duplicated'] == 0]
    filteredDf.drop(columns=['duplicated'],inplace=True)
    
    return(filteredDf)

def main(epiFile):

#   input_file = f'{epiPath}/trainingCombinedEpi.epi.cc.summary'
    input_file = epiFile
    
    output_file = f'{epiFile}.filtered'
    
    # Create a set to track unique pairs
    seen_pairs = set()
    lines_to_keep = []
    
    # Read the original file to preserve exact formatting
    with open(input_file, 'r') as f:
        header = f.readline()  # Save header
        lines_to_keep.append(header)
        
        for line in f:
            fields = line.strip().split()
            snp1, snp2 = fields[1], fields[-1]  # Assuming SNP and BEST_SNP are the 1st and last fields
            
            # Create a sorted pair as a key
            pair_key = tuple(sorted([snp1, snp2]))
            
            # If we haven't seen this pair yet, keep the line
            if pair_key not in seen_pairs:
                seen_pairs.add(pair_key)
                lines_to_keep.append(line)
                
    # Write the filtered lines to output file
    with open(output_file, 'w') as f:
        f.writelines(lines_to_keep)
        
    print(f"Filtered file saved as {output_file}")

    
    
if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description="removing redundant epi pairs that are reversed...")
    parser.add_argument("--epi_file", help="Path to the input epi folder")

    
    args = parser.parse_args()

    # Prefer command-line input if provided; fallback to env var
#   epi_file = '/Users/kerimulterer/prsInteractive/results/celiacDisease/epiFiles/trainingCombinedEpi.epi.cc.summary'
    epi_file = args.epi_file or os.environ.get("EPI_FILE")
    print(f"[PYTHON] Reading from: {epi_file}")
    
    
    if not epi_file:
        raise ValueError("You must provide epi file path via --epi_file or set the EPI_FILE environment variable.")

        
    main(epi_file)
    