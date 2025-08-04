#!/usr/bin/env python3

import pandas as pd
import os
import argparse

def add_iid_simple(long_df, phenotype_df):
    
    """
    Simplest approach - directly repeat phenotype data to match long_df length.
    """
    
    n_models = len(long_df['model'].unique())
    n_individuals = len(phenotype_df)
    
    print(f"Adding IID for {n_models} models × {n_individuals} individuals")
    
    # Create result DataFrame
    result_df = long_df.copy()
    
    # Simply repeat the phenotype data n_models times
    result_df['IID'] = list(phenotype_df['IID']) * n_models

    
    return result_df



def add_row_efficient(output_path, new_df):
    """
    Efficiently add row to CSV file without reading entire file.
    """
    
    if os.path.exists(output_path):
        # File exists - append without reading full file
        print(f"Appending to existing file: {output_path}")
        
        # Append to file (no header since file exists)
        new_df.to_csv(output_path, mode='a', header=False, index=False)
        
    else:
        # File doesn't exist - create new
        print(f"Creating new file: {output_path}")
            
        # Save with header
        new_df.to_csv(output_path, index=False)

def calculate_precision_recall(fp, tp, fn):
    if (fp + tp) == 0:
        precision = 0
    else:
        precision = tp / (fp + tp)
        
    if (fn + tp) == 0:
        recall = 0
    else:
        recall = tp / (fn + tp)
        
    return (precision, recall)

def calculate_fp_fn_tp_percentile(df, cohort, phenotype_col='PHENOTYPE'):
    """
    Calculate metrics for top 20% (high risk) vs bottom 80% (low risk)
    Assumes phenotype: 1 = case, 0 = control
    """
    # Calculate 80th percentile threshold
    threshold_80 = df[cohort].quantile(0.8)
    
    # Define high risk (top 20%) and low risk (bottom 80%)
    hr = df[df[cohort] >= threshold_80]  # Top 20% (high risk)
    lr = df[df[cohort] < threshold_80]   # Bottom 80% (low risk)
    
    print(f"Threshold (80th percentile): {threshold_80:.4f}")
    print(f"High risk group size: {len(hr)} ({len(hr)/len(df)*100:.1f}%)")
    print(f"Low risk group size: {len(lr)} ({len(lr)/len(df)*100:.1f}%)")
    
    # Calculate confusion matrix components
    # Assuming: 1 = case, 0 = control
    tp = hr[hr[phenotype_col] == 2].shape[0]  # True positives: cases in high risk
    fp = hr[hr[phenotype_col] == 1].shape[0]  # False positives: controls in high risk
    fn = lr[lr[phenotype_col] == 2].shape[0]  # False negatives: cases in low risk
    tn = lr[lr[phenotype_col] == 1].shape[0]  # True negatives: controls in low risk
    
    precision, recall = calculate_precision_recall(fp, tp, fn)
    
    # Additional metrics
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    sensitivity = recall  # Same as recall
    
    print(f"TP: {tp}, FP: {fp}, FN: {fn}, TN: {tn}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall/Sensitivity: {recall:.4f}")
    print(f"Specificity: {specificity:.4f}")
    
    return {'false_pos':fp, 'false_neg': fn, 'true_pos': tp, 'true_neg': tn, 'precision':precision, 'recall':recall, 'specificity': specificity, 'sensitivity': sensitivity}

# Alternative version if your phenotypes are coded as 1=control, 2=case
def calculate_fp_fn_tp_percentile_alt_coding(df, cohort, phenotype_col='PHENOTYPE'):
    """
    Version for phenotype coding: 1 = control, 2 = case
    """
    threshold_80 = df[cohort].quantile(0.8)
    
    hr = df[df[cohort] >= threshold_80]  # Top 20%
    lr = df[df[cohort] < threshold_80]   # Bottom 80%
    
    # For 1=control, 2=case coding
    tp = hr[hr[phenotype_col] == 2].shape[0]  # Cases in high risk
    fp = hr[hr[phenotype_col] == 1].shape[0]  # Controls in high risk  
    fn = lr[lr[phenotype_col] == 2].shape[0]  # Cases in low risk
    tn = lr[lr[phenotype_col] == 1].shape[0]  # Controls in low risk
    
    precision, recall = calculate_precision_recall(fp, tp, fn)
    
    # Additional metrics
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    sensitivity = recall  # Same as recall
    
    return {'false_pos':fp, 'false_neg': fn, 'true_pos': tp, 'true_neg': tn, 'precision':precision, 'recall':recall, 'specificity': specificity, 'sensitivity': sensitivity}
        
def calculate_cases_exclusive(df, cohort, threshold_percentile=80, case_value=1,main_col='main'):
    """
    Calculate cases found exclusively by a cohort compared to others.
    
    Args:
        df: DataFrame with PRS scores and PHENOTYPE
        cohort: Column name of the cohort to analyze
        threshold_percentile: Percentile threshold for high-risk (default 80th = top 20%)
        case_value: Value representing cases in PHENOTYPE column (usually 1)
    """
    
    print(f"\n=== Analyzing exclusive cases for {cohort} ===")
    
    # Calculate dynamic thresholds instead of hardcoded 8
    main_threshold = df[main_col].quantile(threshold_percentile / 100)
    cohort_threshold = df[cohort].quantile(threshold_percentile / 100)
    
    print(f"Main threshold (top {100-threshold_percentile}%): {main_threshold:.4f}")
    print(f"{cohort} threshold (top {100-threshold_percentile}%): {cohort_threshold:.4f}")
    
    # Count cases found with main model
    n_main = df[(df[main_col] >= main_threshold) & (df['PHENOTYPE'] == case_value)].shape[0]
    print(f"Cases found with main model: {n_main}")
    
    # Count high-risk cases for current cohort
    hr_cohort = df[(df[cohort] >= cohort_threshold) & (df['PHENOTYPE'] == case_value)]
    print(f"High-risk cases with {cohort}: {hr_cohort.shape[0]}")
    
    # Count cases found by cohort but missed by main
    extra_cases_vs_main = df[
        (df[main_col] < main_threshold) & 
        (df[cohort] >= cohort_threshold) & 
        (df['PHENOTYPE'] == case_value)
    ].shape[0]
    
    print(f"Extra cases found by {cohort} vs main: {extra_cases_vs_main}")
    
    # Get all PRS columns (excluding PHENOTYPE and other non-PRS columns)
    prs_columns = [col for col in df.columns 
                   if col not in ['PHENOTYPE', 'IID'] 
                   and not col.startswith('PRScr')  # Exclude any score columns
                   and col != cohort]  # Exclude current cohort
    
    print(f"Comparing against PRS models: {prs_columns}")
    
    # Find cases that are high-risk ONLY in the current cohort
    # (low-risk in ALL other PRS models)
    condition_low_in_others = pd.Series([True] * len(df))
    
    for other_col in prs_columns:
        if other_col in df.columns:
            other_threshold = df[other_col].quantile(threshold_percentile / 100)
            condition_low_in_others = condition_low_in_others & (df[other_col] < other_threshold)
            
    # Cases found exclusively by current cohort
    unique_cases = df[
        condition_low_in_others & 
        (df[cohort] >= cohort_threshold) & 
        (df['PHENOTYPE'] == case_value)
    ]
    
    n_unique_cases = unique_cases.shape[0]
    
    # Calculate improvement percentages
    percent_improvement_vs_main = calculate_percent_improvement(extra_cases_vs_main, n_main)
    percent_improvement_unique = calculate_percent_improvement(n_unique_cases, n_main)
    
    print(f"Cases found ONLY by {cohort} (missed by all others): {n_unique_cases}")
    print(f"Improvement vs main: {percent_improvement_vs_main:.2f}%")
    print(f"Unique contribution: {percent_improvement_unique:.2f}%")
    
    return {
        'cohort': cohort,
        'main_cases': n_main,
        'cohort_high_risk_cases': hr_cohort.shape[0],
        'extra_vs_main': extra_cases_vs_main,
        'unique_cases': n_unique_cases,
        'improvement_vs_main_pct': percent_improvement_vs_main,
        'unique_contribution_pct': percent_improvement_unique,
        'main_threshold': main_threshold,
        'cohort_threshold': cohort_threshold
    }
    
def calculate_percent_improvement(extra_cases, baseline_cases):
    """Calculate percentage improvement."""
    if baseline_cases == 0:
        return 0 if extra_cases == 0 else float('inf')
    return (extra_cases / baseline_cases) * 100


def calculate_precision_recall_improvement(scoresPath,model_type='prs'):
    outputPath = f'{scoresPath}/model_recall_precision_improvement.csv'

    if model_type == "prs":
        for h in ['validation','holdout',]:
            if h == 'validation':
                filePath = f'{scoresPath}/combinedPRSGroups.csv'
            else:
                filePath = f'{scoresPath}/CombinedPRSGroups.holdout.csv'
            df = pd.read_csv(filePath)
            cohorts = [col for col in df.columns if 'scaled_prs' in col]

            for cohort in cohorts:
                fp_fn_dict = calculate_fp_fn_tp_percentile_alt_coding(df[cohorts+['PHENOTYPE']],cohort)
                performance_values = calculate_cases_exclusive(df, cohort, threshold_percentile=80, case_value=2,main_col='scaled_prs_main')
                tempDf = pd.DataFrame(fp_fn_dict | performance_values,index=[0])
                tempDf['data_type'] = 'prs'
                tempDf['data_subset'] = h
                add_row_efficient(outputPath, tempDf)
    else: #calculate stats for trained models
        filePath = f'{scoresPath}/predictProbsReducedFinalModel.csv'
        phenotype = pd.read_csv(f'{scorePath}/combinedPRSGroups.csv',usecols=['IID','PHENOTYPE'])
        df = pd.read_csv(filePath)
        
        #remove covariate model
        df = df[df['model'] != 'covariate']
        
        if 'IID' not in df.columns:
            df = add_iid_simple(df, phenotype)

        # Convert to wide format
        df_wide = df.pivot(index=['IID'], columns='model', values='yProba')
        df_wide.reset_index(inplace=True)
        
        # Flatten column names (remove multi-level index)
        df_wide.columns.name = None
        
        df_wide = df_wide.merge(phenotype,on='IID',how='left')
        
        #cohorts in the model column
        cohorts = df['model'].unique()
        

        
        for cohort in cohorts:
            fp_fn_dict = calculate_fp_fn_tp_percentile_alt_coding(df_wide[cohorts.tolist()+['PHENOTYPE']],cohort)
#           missed_list = calculate_cases_exclusive(df_wide[cohorts.tolist()+['PHENOTYPE']],cohort)
            performance_values = calculate_cases_exclusive(df_wide, cohort, threshold_percentile=80, case_value=2)
            tempDf = pd.DataFrame(fp_fn_dict | performance_values,index=[0])
            tempDf['data_type'] = 'penalized_model'
            tempDf['data_subset'] = 'validation'
            add_row_efficient(outputPath, tempDf)
        

            

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="calculating model performance of G vs other using predictions of trained models and PRS calculations....")
    parser.add_argument("--pheno_path", help="Path to the input pheno folder")


    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_path or set the PHENO_PATH environment variable.")
    

    scorePath = f'{pheno_path}/scores'
    for t in ['model','prs']:
        calculate_precision_recall_improvement(scorePath,model_type=t)