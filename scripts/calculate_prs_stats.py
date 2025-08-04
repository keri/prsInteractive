#!/usr/bin/env python3

import pandas as pd
import os





def add_row_efficient(output_path, new_row_data, columns=None):
    """
    Efficiently add row to CSV file without reading entire file.
    """
    
    if os.path.exists(output_path):
        # File exists - append without reading full file
        print(f"Appending to existing file: {output_path}")
        
        # Create DataFrame for new row
        if isinstance(new_row_data, dict):
            new_df = pd.DataFrame([new_row_data])
        else:
            new_df = pd.DataFrame([new_row_data], columns=columns)
            
        # Append to file (no header since file exists)
        new_df.to_csv(output_path, mode='a', header=False, index=False)
        
    else:
        # File doesn't exist - create new
        print(f"Creating new file: {output_path}")
        
        # Create DataFrame
        if isinstance(new_row_data, dict):
            df = pd.DataFrame([new_row_data])
        else:
            df = pd.DataFrame([new_row_data], columns=columns)
            
        # Save with header
        df.to_csv(output_path, index=False)

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
    
    return [fp, fn, tp, tn, precision, recall, specificity, sensitivity]

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
    
    return [fp, fn, tp, tn, precision, recall]
        
def calculate_percent_improvement(extra_cases,n_main):
    return(extra_cases/(n_main+extra_cases))
    
def calculate_cases_exclusive(df,cohorts,cohort):
    main_col = [col for col in df.columns if '_main' in col][0]
    
    #n cases found with main
    n_main = df[(df[main_col] > 8) & (df['PHENOTYPE'] == 2)].shape[0]
    
    #get the hr cases for cohort
    hr = df[(df[cohort] > 8) & (df['PHENOTYPE'] == 2)]
    #count number of these not found in main
    extra_cases = df[(df[main_col] < 9) & (df[cohort] > 8) & (df['PHENOTYPE'] == 2)].shape[0]
    
    #percent improvement compared to main
    percent_improvement = calculate_percent_improvement(extra_cases, n_main)
    
    print('number of high risk cases = ',hr.shape[0])
    temp = hr.drop(columns=[cohort])
    print(f'number missed with main and found with {cohort} = {extra_cases}')
    rest_of_cohorts = [col for col in temp.columns if col != 'PHENOTYPE']
    rest_of_cohorts = [col for col in rest_of_cohorts if 'PRScr' not in col]
    for c in rest_of_cohorts:
        temp = temp[temp[c] < 9]
    
    n_missed_all_others = temp.shape[0]
    print(f'number missed with all others and found with {cohort} = {n_missed_all_others}')
        
    return([n_missed_all_others,extra_cases,percent_improvement])
        
def calculate_precision_recall_improvement(scorePath,model_type='prs'):
    outputPath = f'{scoresPath}/model_recall_precision_improvement.csv'
    stats_columns = ['model_type','model','data','false_positive','false_negative','true_positive','true_negative','precision','recall','n_missed_all_others','n_missed_with_main','percent_improvement_to_main']
    if model_type == "prs":
        for h in ['test','holdout',]:
            if h == 'test':
                filePath = f'{scorePath}/combinedPRSGroups.csv'
            else:
                filePath = f'{scorePath}/CombinedPRSGroups.holdout.csv'
            df = pd.read_csv(filePath)
            cohorts = [col for col in df.columns if 'scaled_prs' in col]

            for cohort in cohorts:
                fp_fn_list = calculate_fp_fn_tp(df[cohorts+['PHENOTYPE']],cohort)
                missed_list = calculate_cases_exclusive(df[cohorts+['PHENOTYPE']],cohorts,cohort)
                values = fp_fn_list + missed_list
                row = [model_type,cohort.split('_')[-1],h]+values
                add_row_efficient(outputPath, row, columns=stats_columns)
    else: #calculate stats for trained models
        filePath = f'{scoresPath}/predictProbsReducedFinalModel.csv'
        phenotype = pd.read_csv(f'{scorePath}/combinedPRSGroups.csv',usecols=['IID','PHENOTYPE'])
        df = pd.read_csv(filePath)
        
        #cohorts in the model column
        cohorts = df['model'].unique()
        
        df = df.merge(phenotype,on='IID',how='left')
        
        # Convert to wide format
        df_wide = df.pivot(index=['IID', 'PHENOTYPE'], columns='model', values='yProba')
    
        # Reset index to make IID and PHENOTYPE regular columns
        df_wide = df_wide.reset_index()
        
        # Flatten column names (remove multi-level index)
        df_wide.columns.name = None
        
        for cohort in cohorts:
            fp_fn_list = calculate_fp_fn_tp(df[cohorts+['PHENOTYPE']],cohort)
            missed_list = calculate_cases_exclusive(df[cohorts+['PHENOTYPE']],cohorts,cohort)
            values = fp_fn_list + missed_list
            row = [model_type,cohort.split('_')[-1],'test']+values
            add_row_efficient(outputPath, row, columns=stats_columns)
        

            

if __name__ == '__main__':
    
    pheno = 'celiacDisease'
    scorePath = f'/Users/kerimulterer/prsInteractive/results/{pheno}/scores'
    calculate_precision_recall_improvement(scorePath,type='model')