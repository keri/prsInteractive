#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import sys

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import sys

def main(df2, pheno, pheno_path, icd_code="E11", pheno_str="type 2 diabetes"):
    '''
    input:
        df2 : pandas.DataFrame() with 
        minimum columns needed = ['Non-cancer illness code, self-reported | Instance 0',
        'Non-cancer illness code, self-reported | Instance 1',
        'Non-cancer illness code, self-reported | Instance 2',
        'Non-cancer illness code, self-reported | Instance 3',
        'Diagnoses - main ICD10',
        'Participant ID'

    output:
        Four .txt files space delimited
             trainingID.txt, testID.txt, holdoutID.txt, pheno.txt
    '''
    
    print(f"Processing phenotype: {pheno}")
    print(f"Using ICD code: {icd_code}")
    print(f"Using phenotype string: {pheno_str}")
    print(f"Output path: {pheno_path}")
    print(f"Input dataframe shape: {df2.shape}")
    
    # Validate inputs
    if not icd_code or not isinstance(icd_code, str):
        print(f"ERROR: Invalid ICD code: {icd_code}. Must be a non-empty string.")
        return False
    
    if not pheno_str or not isinstance(pheno_str, str):
        print(f"ERROR: Invalid phenotype string: {pheno_str}. Must be a non-empty string.")
        return False
    
    # Print column names for debugging
    print("Available columns:")
    for i, col in enumerate(df2.columns):
        print(f"  {i}: {col}")
        
    # Initialize phenotype column with controls (1)
    df2['phenotype'] = 1
    
    ###########################  DEFINE PHENOTYPE  #######################
    try:
        if pheno == "type2Diabetes":
            # Create boolean mask for cases
            mask = (
                (df2['Date E11 first reported (non-insulin-dependent diabetes mellitus)'].notna() & 
                 (df2['Date E11 first reported (non-insulin-dependent diabetes mellitus)'] != '')) |
                (df2['Source of report of E11 (non-insulin-dependent diabetes mellitus)'].notna() & 
                 (df2['Source of report of E11 (non-insulin-dependent diabetes mellitus)'] != '')) |
                (df2['Non-cancer illness code, self-reported | Instance 0'].str.contains('type 2 diabetes', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 1'].str.contains('type 2 diabetes', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 2'].str.contains('type 2 diabetes', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 3'].str.contains('type 2 diabetes', na=False)) |
                (df2['Diagnoses - main ICD10'].str.contains('E11', na=False))
            )
        elif pheno == 'celiacDisease':
            mask = (
                (df2['Date K90 first reported (intestinal malabsorption)'].notna() & 
                 (df2['Date K90 first reported (intestinal malabsorption)'] != '')) |
                (df2['Source of report of K90 (intestinal malabsorption)'].notna() & 
                 (df2['Source of report of K90 (intestinal malabsorption)'] != '')) |
                (df2['Non-cancer illness code, self-reported | Instance 0'].str.contains('coeliac disease', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 1'].str.contains('coeliac disease', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 2'].str.contains('coeliac disease', na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 3'].str.contains('coeliac disease', na=False)) |
                (df2['Diagnoses - main ICD10'].str.contains('K90', na=False))
            )
        else:
            # Generic phenotype matching - ensure strings are valid
            mask = (
                (df2['Non-cancer illness code, self-reported | Instance 0'].str.contains(pheno_str, na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 1'].str.contains(pheno_str, na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 2'].str.contains(pheno_str, na=False)) |
                (df2['Non-cancer illness code, self-reported | Instance 3'].str.contains(pheno_str, na=False)) |
                (df2['Diagnoses - main ICD10'].str.contains(icd_code, na=False))
            )
            
        # Set cases to 2
        df2.loc[mask, 'phenotype'] = 2
        
        # Count cases and controls
        case_count = (df2['phenotype'] == 2).sum()
        control_count = (df2['phenotype'] == 1).sum()
        
        print(f'Total cases found: {case_count}')
        print(f'Total controls: {control_count}')
        print(f'Phenotype prevalence: {case_count/df2.shape[0]:.4f}')
        
        if case_count == 0:
            print("ERROR: No cases found! Check phenotype definition.")
            print("Sample of relevant columns:")
            relevant_cols = [col for col in df2.columns if any(term in col.lower() for term in ['illness', 'icd10', 'diagnos'])]
            if relevant_cols:
                print(df2[relevant_cols[:5]].head())
            return False
        
    except KeyError as e:
        print(f"ERROR: Missing column {e}")
        print("Available columns that might be relevant:")
        relevant_cols = [col for col in df2.columns if any(term in col.lower() for term in ['illness', 'icd', 'diagnos', 'diabetes'])]
        for col in relevant_cols:
            print(f"  - {col}")
        return False
    
    ##################### SPLIT DATA INTO TRAIN/TEST/HOLDOUT  #################################
    
    # Ensure we have both cases and controls for splitting
    cases = df2[df2['phenotype'] == 2]
    controls = df2[df2['phenotype'] == 1]
    
    if len(cases) < 10:
        print(f"WARNING: Very few cases ({len(cases)}). Consider adjusting phenotype definition.")
        
    # Split data maintaining class proportions
    training_cases = cases.sample(frac=0.70, random_state=1)
    training_controls = controls.sample(frac=0.70, random_state=1)
    training = pd.concat([training_cases, training_controls])
    
    print(f'Training data consists of {training.shape[0]} people')
    
    # Get remaining data
    remaining_cases = cases[~cases.index.isin(training_cases.index)]
    remaining_controls = controls[~controls.index.isin(training_controls.index)]
    
    # Split remaining data between test and holdout (roughly 2:1 ratio)
    test_cases = remaining_cases.sample(frac=0.6666, random_state=1)
    test_controls = remaining_controls.sample(frac=0.6666, random_state=1)
    test = pd.concat([test_cases, test_controls])
    
    print(f'Test data consists of {test.shape[0]} people')
    
    # Holdout is the remainder
    holdout_cases = remaining_cases[~remaining_cases.index.isin(test_cases.index)]
    holdout_controls = remaining_controls[~remaining_controls.index.isin(test_controls.index)]
    holdout = pd.concat([holdout_cases, holdout_controls])
    
    print(f'Holdout data consists of {holdout.shape[0]} people')
    
    # Print statistics for each split
    for split_name, split_data in [('training', training), ('test', test), ('holdout', holdout)]:
        cases_in_split = (split_data['phenotype'] == 2).sum()
        controls_in_split = (split_data['phenotype'] == 1).sum()
        case_percentage = (cases_in_split / len(split_data)) * 100
        
        print(f'\n{split_name.capitalize()} set statistics:')
        print(f'  Cases: {cases_in_split}')
        print(f'  Controls: {controls_in_split}')
        print(f'  Case percentage: {case_percentage:.2f}%')
        
    # Create output DataFrames
    trainingID = training[['Participant ID']].copy()
    trainingID['IID'] = training['Participant ID']
    
    holdoutID = holdout[['Participant ID']].copy()
    holdoutID['IID'] = holdout['Participant ID']
    
    testID = test[['Participant ID']].copy()
    testID['IID'] = test['Participant ID']
    
    phenotype = df2[['Participant ID']].copy()
    phenotype['IID'] = df2['Participant ID']
    phenotype['phenotype'] = df2['phenotype']
    
    # Ensure output directory exists
    os.makedirs(pheno_path, exist_ok=True)
    
    # Save files
    try:
        testID.to_csv(f'{pheno_path}/testID.txt', sep=' ', header=False, index=False)
        holdoutID.to_csv(f'{pheno_path}/holdoutID.txt', sep=' ', header=False, index=False)
        trainingID.to_csv(f'{pheno_path}/trainingID.txt', sep=' ', header=False, index=False)
        phenotype.to_csv(f'{pheno_path}/pheno.txt', sep=' ', header=False, index=False)
        
        print(f"\nSuccessfully saved files to {pheno_path}:")
        for filename in ['testID.txt', 'holdoutID.txt', 'trainingID.txt', 'pheno.txt']:
            filepath = os.path.join(pheno_path, filename)
            if os.path.exists(filepath):
                print(f"  ✅ {filename} ({os.path.getsize(filepath)} bytes)")
            else:
                print(f"  ❌ {filename} (not found)")
                
        return True
    
    except Exception as e:
        print(f"ERROR saving files: {e}")
        return False
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Creating phenotype file and train/test/holdout splits")
    parser.add_argument("--data_path", required=False, help="Path to the prsInteractive data folder")
    parser.add_argument("--pheno_path", required=False, help="Path to the output phenotype results folder")
    parser.add_argument("--pheno", required=False, help="Phenotype to analyze")
    parser.add_argument("--pheno_str", required=False, help="Phenotype string used in Non-cancer illness code, self-reported field")
    parser.add_argument("--icd10", required=False, help="ICD 10 code of phenotype")
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var with defaults
    data_path = args.data_path or os.environ.get("DATA_PATH")
    pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
    pheno = args.pheno or os.environ.get("PHENO", "type2Diabetes")  # Default phenotype
    icd10 = args.icd10 or os.environ.get("ICD10", "E11")  # Default ICD code
    pheno_str = args.pheno_str or os.environ.get("PHENO_STR", "type 2 diabetes")  # Default string
    
    print(f"[PYTHON] Reading from: {data_path}")
    print(f"[PYTHON] Output to: {pheno_path}")
    print(f"[PYTHON] Phenotype: {pheno}")
    print(f"[PYTHON] ICD code: {icd10}")
    print(f"[PYTHON] Phenotype string: {pheno_str}")
    
    # Validate required parameters
    if not data_path:
        print("ERROR: DATA_PATH not provided via command line or environment variable")
        sys.exit(1)
        
    if not pheno_path:
        print("ERROR: PHENO_PATH not provided via command line or environment variable")
        sys.exit(1)
        
    # Check if data path exists
    if not os.path.exists(data_path):
        print(f"ERROR: Data path {data_path} does not exist!")
        sys.exit(1)
        
    # Check if participant.csv exists
    participant_file = os.path.join(data_path, 'participant.csv')
    if not os.path.exists(participant_file):
        print(f"ERROR: participant.csv not found at {participant_file}")
        print(f"Available files in {data_path}:")
        try:
            for f in os.listdir(data_path):
                print(f"  - {f}")
        except:
            print("  (cannot list directory)")
        sys.exit(1)
        
    print(f"Loading data from {participant_file}...")
    
    try:
        # Read participant dataset
        df = pd.read_csv(participant_file)
        print(f"Loaded dataset with shape: {df.shape}")
        
        # Fill NaN values with empty strings for string operations
        df2 = df.fillna('')
        
        # Run main processing
        success = main(df2, pheno, pheno_path, icd_code=icd10, pheno_str=pheno_str)
        
        if success:
            print("\n✅ Phenotype processing completed successfully!")
            sys.exit(0)
        else:
            print("\n❌ Phenotype processing failed!")
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)