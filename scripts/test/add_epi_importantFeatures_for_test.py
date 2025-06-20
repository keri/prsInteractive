#!/usr/bin/env python3

import pandas as pd
import os
import argparse



def main(pheno_path):
    # add a few epi features to test the workflow
    
    df = pd.DataFrame()
    df['feature'] = ['disease_66,disease_88','disease_88,disease_66','nullB_78,disease_68','disease_24,disease_29']
    df['shap_zscore'] = [2.5,2.5,3.0,7.5]
    df['data_type'] = 'epi'
    
    with open(f'{pheno_path}/scores/importantFeaturesPostShap.csv',mode='a',newline='') as f:
        df.to_csv(f,index=False, header=False)
        f.close()
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="creating LD SNP list ...")
    parser.add_argument("--pheno_folder", help="Path to the input pheno")
    
    args = parser.parse_args()
    
    pheno_path = args.pheno_folder or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    #pheno_path = '/Users/kerimulterer/prsInteractive/results/type2Diabetes_test'
    
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_folder or set the PHENO_PATH environment variable.")
        

    
    main(pheno_path)