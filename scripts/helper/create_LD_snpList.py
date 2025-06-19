#!/usr/bin/env python3

import pandas as pd
import os
import argparse




# add a few epi features to test the workflow

def	main(phenoPath):
    
    df = pd.read_csv(f'{phenoPath}/scores/importantFeaturesPostShap.csv')
    
    mainFeatures = df[df['data_type'] == 'main']['feature'].tolist()
    
    epiFeatures = df[df['data_type'] == 'epi']['feature'].tolist()
    
    epiFeatures2 = []
    for feature in epiFeatures:
        featureList = feature.split(',')
        for f in featureList:
            epiFeatures2.append(f)
    LDSnps = pd.DataFrame({'feature':list(set(epiFeatures2+mainFeatures))})
    LDSnps.to_csv(f'{phenoPath}/finalModelLDSnps.txt',header=False,index=False)

    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="creating LD SNP list ...")
    parser.add_argument("--pheno_path", help="Path to the input pheno")
    
    args = parser.parse_args()
    
    pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
#   pheno_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes'
    
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_path or set the PHENO_PATH environment variable.")
    

    
    main(pheno_path)