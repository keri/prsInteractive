#!/usr/bin/env python3

import pandas as pd
import os
import argparse




# add a few epi features to test the workflow

def	main(phenoPath,feature_scores_file):
    
    df = pd.read_csv(feature_scores_file)
    
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
    parser.add_argument("--pheno_data", help="Path to the input scores")
    
    
    args = parser.parse_args()
    
    pheno_data = args.pheno_path or os.environ.get("PHENO_DATA")
#   print(f"[PYTHON] Reading from: {pheno_path}")

    
#   pheno_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes'
    feature_scores_file = f"{pheno_data}/scores/importantFeaturesPostShap.csv"
    print(f"[PYTHON] Reading features for LD from: {feature_scores_file}")
    
    if not pheno_data:
        raise ValueError("You must provide a data scores path via --pheno_data or set the PHENO_DATA environment variable.")
        
    
    main(pheno_data,feature_scores_file)

    

    