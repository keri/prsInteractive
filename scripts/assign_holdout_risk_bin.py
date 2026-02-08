#!/usr/bin/env python3

import pandas as pd
import numpy as np

from helper.data_wrangling import *
from helper.calculate_prs import *
    
def create_case_control_histogram(df,pheno_col,continous_col,figPath,figsize=(12,6)):
    """
    Plot the distribution of a cases and controls
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing prs/bin calculations
    pheno_col : str
        Column name of the binarized version of the phenotype
    continous_col : str
        column name of the measure to be plotted
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        Distribution plot
    """
    plt.figure(figsize=figsize)
    
    if not set(df[pheno_col].unique()).issubset({0, 1}):
        df['phenotype'] = df['PHENOTYPE'] - 1
        pheno_col = 'phenotype'
        
        
    # Split data
    cases = df[df[pheno_col] == 1][[continous_col]]
    controls = df[df[pheno_col] == 0][[continous_col]]
    
    # Plot histograms
    plt.hist(controls, bins=30, alpha=0.6, label='Controls', color='skyblue', edgecolor='black')
    plt.hist(cases, bins=30, alpha=0.6, label='Cases', color='salmon', edgecolor='black')
    
    case_mean = cases[continous_col].mean()
    # Plot threshold line
    plt.axvline(x=case_mean, color='black', linestyle='--', linewidth=2)
    
    
    plt.title(f'Distribution of cases/controls across: {continous_col}')
    plt.xlabel(continous_col)
    plt.ylabel('Frequency')
    plt.grid(alpha=0.3)
    
    plt.savefig(f'{figPath}/histogramCaseControl.{continous_col}.png')
    plt.close()
    


def create_risk_bins(training_df, n_bins=1000, phenotype_col='PHENOTYPE',use_epi_main=False):
    """
    Create risk bins for PRS scores and calculate bin statistics for cases and controls separately.
    
    Parameters:
    training_df: DataFrame with PRS columns, phenotype column, and IID as index
    n_bins: Number of risk bins to create (default 1000)
    phenotype_col: Name of the column containing case/control status (default 'status')
                   Assumes 1 = case, 0 = control
    
    Returns:
    bin_stats: Dictionary containing bin statistics for each PRS type
    """
    if 2 in training_df[phenotype_col].unique():
        training_df['phenotype'] = training_df['PHENOTYPE'] - 1
        phenotype_col = 'phenotype'
        
    if use_epi_main:
        prs_columns = ['scaled_prs_main', 'scaled_prs_epi','scaled_prs_cardio','scaled_prs_epi+main']
    else:
        prs_columns = ['scaled_prs_main', 'scaled_prs_epi','scaled_prs_cardio']
    
    bin_stats = {}
    
    print(f"Creating {n_bins} risk bins for each PRS type...")
    print(f"Training data shape: {training_df.shape}")
    
    # Check if phenotype column exists
    if phenotype_col not in training_df.columns:
        raise ValueError(f"Phenotype column '{phenotype_col}' not found in training data")
        
    # Calculate overall case rate for reference
    overall_case_rate = training_df[phenotype_col].mean()
    print(f"Overall case rate: {overall_case_rate:.4f}")
    
    for prs_col in prs_columns:
        if prs_col not in training_df.columns:
            print(f"Warning: {prs_col} not found in training data. Skipping...")
            continue
        
        print(f"\nProcessing {prs_col}...")
        
        # Create bins using quantiles to ensure equal sample sizes
        # Use only rows with non-missing PRS values
        valid_mask = training_df[prs_col].notna()
        prs_values = training_df.loc[valid_mask, prs_col]
        phenotype_values = training_df.loc[valid_mask, phenotype_col]
        
        # Create bin edges using quantiles
        bin_edges = np.linspace(0, 1, n_bins + 1)
        quantile_edges = prs_values.quantile(bin_edges).values
        
        # Assign bin numbers (0 to n_bins-1)
        bin_assignments = pd.cut(prs_values, bins=quantile_edges, 
                               labels=range(n_bins), include_lowest=True)
        
        # Calculate bin statistics for cases and controls separately
        bin_data = []
        for bin_num in range(1,n_bins+1):
            bin_mask = bin_assignments == bin_num
            bin_prs_values = prs_values[bin_mask]
            bin_phenotype_values = phenotype_values[bin_mask]
            
            if len(bin_prs_values) > 0:
                # Separate cases and controls
                cases_mask = bin_phenotype_values == 1
                controls_mask = bin_phenotype_values == 0
                
                n_cases = cases_mask.sum()
                n_controls = controls_mask.sum()
                n_total = len(bin_prs_values)
                
                # Calculate case rate in this bin
                case_rate = n_cases / n_total if n_total > 0 else np.nan
                
                # Calculate odds ratio relative to overall case rate
                if overall_case_rate > 0 and overall_case_rate < 1 and case_rate > 0 and case_rate < 1:
                    odds_ratio = (case_rate / (1 - case_rate)) / (overall_case_rate / (1 - overall_case_rate))
                else:
                    odds_ratio = np.nan
                    
                bin_stats_dict = {
                    'bin_number': bin_num,
                    'min_value': bin_prs_values.min(),
                    'max_value': bin_prs_values.max(),
                    'mean_value': bin_prs_values.mean(),
                    'n_cases': int(n_cases),
                    'n_controls': int(n_controls),
                    'n_total': n_total,
                    'case_rate': case_rate,
                    'odds_ratio': odds_ratio,
                    'percentile_min': (bin_num / n_bins) * 100,
                    'percentile_max': ((bin_num + 1) / n_bins) * 100
                }
            else:
                # Handle empty bins
                bin_stats_dict = {
                    'bin_number': bin_num,
                    'min_value': np.nan,
                    'max_value': np.nan,
                    'mean_value': np.nan,
                    'n_cases': 0,
                    'n_controls': 0,
                    'n_total': 0,
                    'case_rate': np.nan,
                    'odds_ratio': np.nan,
                    'percentile_min': (bin_num / n_bins) * 100,
                    'percentile_max': ((bin_num + 1) / n_bins) * 100
                }
                
            bin_data.append(bin_stats_dict)
            
        # Convert to DataFrame for easier handling
        bin_df = pd.DataFrame(bin_data)
        
        # Store statistics
        bin_stats[prs_col] = {
            'bin_dataframe': bin_df,
            'bin_edges': quantile_edges,
            'total_bins': n_bins,
            'overall_case_rate': overall_case_rate
        }
        
        print(f"  Created {n_bins} bins")
        print(f"  Value range: {prs_values.min():.4f} to {prs_values.max():.4f}")
        print(f"  Non-empty bins: {bin_df['n_total'].gt(0).sum()}")
        print(f"  Total cases: {bin_df['n_cases'].sum()}")
        print(f"  Total controls: {bin_df['n_controls'].sum()}")
        
    return bin_stats


def assign_risk_thresholds(holdout_df, bin_stats, training_scalers=None, threshold_percentile=80):
    """
    Assign risk thresholds to holdout data based on training bin statistics.
    Uses case/control bin characteristics to create a combined risk score.
    
    Parameters:
    holdout_df: DataFrame with PRS columns (may be unscaled)
    bin_stats: Dictionary from create_risk_bins() with case/control statistics
    training_scalers: Dictionary of fitted StandardScaler objects (optional)
    threshold_percentile: Percentile threshold for high-risk classification
    
    Returns:
    holdout_processed: DataFrame with assigned risk bins and thresholds
    """
    
    print(f"\nAssigning risk thresholds to holdout data...")
    print(f"Holdout data shape: {holdout_df.shape}")
    print(f"Using {threshold_percentile}th percentile as high-risk threshold")
    
    holdout_processed = holdout_df.copy()
    
    # Scale holdout data if scalers provided
    if training_scalers is not None:
        print("Scaling holdout data using training scalers...")
        for prs_col in training_scalers.keys():
            if prs_col in holdout_df.columns:
                # Get the base column name (remove 'scaled_' prefix if present)
                base_col = prs_col.replace('scaled_', '') if prs_col.startswith('scaled_') else prs_col
                
                if base_col in holdout_df.columns:
                    scaled_values = training_scalers[prs_col].transform(
                        holdout_df[[base_col]]
                    ).flatten()
                    holdout_processed[f'scaled_{base_col}'] = scaled_values
                    print(f"  Scaled {base_col} -> scaled_{base_col}")
                    
    # Calculate threshold bin (e.g., 800th bin for 80th percentile with 1000 bins)
    threshold_bin = int((threshold_percentile / 100) * 1000)
    print(f"High-risk threshold: bin {threshold_bin} and above")
    
    # Track bins for combined risk calculation
    all_bins = {}  # Store bin assignments for each PRS type
    case_bin_flags = {}  # Track whether each bin is a "case bin"
    
    # Assign bins and thresholds for each PRS type
    for prs_col, stats in bin_stats.items():
        if prs_col not in holdout_processed.columns:
            print(f"Warning: {prs_col} not found in holdout data. Skipping...")
            continue
        
        print(f"\nProcessing {prs_col} in holdout data...")
        
        bin_df = stats['bin_dataframe']
        bin_edges = stats['bin_edges']
        overall_case_rate = stats['overall_case_rate']
        
        # Assign bins to holdout data
        prs_values = holdout_processed[prs_col]
        bin_assignments = pd.cut(prs_values, bins=bin_edges, 
                               labels=range(len(bin_edges)-1), include_lowest=True)
        
        # Convert to numeric (handling NaN)
#       bin_assignments_numeric = pd.to_numeric(bin_assignments, errors='coerce')
        # Convert categorical to integer codes directly (more reliable than to_numeric)
        bin_assignments_numeric = bin_assignments.cat.codes
        
        # Check for unassigned bins (codes = -1 means NaN in categorical)
        n_unassigned = (bin_assignments_numeric == -1).sum()
        if n_unassigned > 0:
            print(f"WARNING: {n_unassigned} values could not be assigned to bins")
            print(f"  PRS range in holdout: [{prs_values.min():.6f}, {prs_values.max():.6f}]")
            print(f"  Bin edges range: [{bin_edges.min():.6f}, {bin_edges.max():.6f}]")
            
            # Assign edge cases to nearest bin
            # Values below minimum edge -> bin 0
            # Values above maximum edge -> last bin
            bin_assignments_numeric = bin_assignments_numeric.copy()
            bin_assignments_numeric[prs_values < bin_edges[0]] = 0
            bin_assignments_numeric[prs_values > bin_edges[-1]] = len(bin_edges) - 2
            
        bin_assignments_numeric = bin_assignments_numeric+1
        
        # Create column names for this PRS type
        base_col = prs_col.replace("scaled_prs_","")
        bin_col = f'{base_col}_centile_bin'
        threshold_col = f'{base_col}_high_risk'
        case_bin_col = f'{base_col}_case_bin'
        
        # Assign bin numbers
        holdout_processed[bin_col] = bin_assignments_numeric
        
        # Determine if each person's bin is a "case bin" (case-enriched)
        # A bin is considered a case bin if its case_rate > overall_case_rate
        is_case_bin = bin_assignments_numeric.map(
            lambda x: bin_df.loc[int(x), 'case_rate'] > overall_case_rate 
            if pd.notna(x) and int(x) < len(bin_df) else False
        )
        holdout_processed[case_bin_col] = is_case_bin.astype(int)
        
        # Assign high-risk status (1 if bin >= threshold_bin, 0 otherwise)
        holdout_processed[threshold_col] = (bin_assignments_numeric >= threshold_bin).astype(int)
        
        # Add threshold value for reference
        threshold_value = bin_df.loc[threshold_bin, 'min_value'] if threshold_bin < len(bin_df) else bin_df['max_value'].max()
        holdout_processed[f'{prs_col}_threshold_value'] = threshold_value
        
        # Store for combined calculation
        all_bins[base_col] = bin_assignments_numeric
        case_bin_flags[base_col] = is_case_bin
        
        # Statistics
        high_risk_count = holdout_processed[threshold_col].sum()
        case_bin_count = is_case_bin.sum()
        total_count = len(holdout_processed.dropna(subset=[prs_col]))
        
        print(f"  Assigned bins: 1-{bin_assignments_numeric.max():.0f}")
        print(f"  High-risk threshold value: {threshold_value:.4f}")
        print(f"  High-risk individuals: {high_risk_count}/{total_count} ({high_risk_count/total_count*100:.1f}%)")
        print(f"  Case-enriched bins: {case_bin_count}/{total_count} ({case_bin_count/total_count*100:.1f}%)")
        
    # Create composite high-risk indicator (high-risk in ANY PRS type)
    risk_columns = [col for col in holdout_processed.columns if col.endswith('_high_risk') and not col.startswith('combined')]
    if risk_columns:
        holdout_processed['any_high_risk'] = holdout_processed[risk_columns].max(axis=1)
        any_high_risk_count = holdout_processed['any_high_risk'].sum()
        print(f"\nComposite high-risk (any PRS type): {any_high_risk_count}/{len(holdout_processed)} ({any_high_risk_count/len(holdout_processed)*100:.1f}%)")
        
    # Calculate combined PRS based on case bin proportion
    print(f"\n{'='*60}")
    print("CALCULATING COMBINED PRS BASED ON CASE BIN PROPORTION")
    print(f"{'='*60}")
    
    if all_bins:
        # Convert to DataFrame for easier calculation
        bins_df = pd.DataFrame(all_bins)
        case_flags_df = pd.DataFrame(case_bin_flags)
        
        # Count number of case bins for each individual
        case_bin_count_per_person = case_flags_df.sum(axis=1)
        
        # Count number of non-NA PRS types for each individual
        valid_prs_count = bins_df.notna().sum(axis=1)
        
        # Calculate proportion of case bins (only for individuals with at least one valid PRS)
        prop_case_bins = case_bin_count_per_person / valid_prs_count.replace(0, np.nan)
        
        # Store case bin metrics
        holdout_processed['n_case_bins'] = case_bin_count_per_person
        holdout_processed['n_valid_prs'] = valid_prs_count
        holdout_processed['prop_case_bins'] = prop_case_bins
        
        # Assign combined bin based on 60% threshold
        combined_bin = np.zeros(len(holdout_processed))
        combined_prs = np.zeros(len(holdout_processed))
        scaled_combined_prs = np.zeros(len(holdout_processed))
        
        # Track assignment reasons for reporting
        n_any_high_risk = 0
        n_prop_high = 0
        n_prop_low = 0
        n_missing = 0
        
        for idx in range(len(holdout_processed)):
            if valid_prs_count.iloc[idx] == 0:
                combined_bin[idx] = np.nan
                combined_prs[idx] = np.nan
                scaled_combined_prs[idx] = np.nan
                n_missing += 1
                
            elif holdout_processed['any_high_risk'].iloc[idx] > 0:
                # Priority 1: If any_high_risk > 0, use MAX bin
                max_bin = bins_df.iloc[idx].max()
                combined_bin[idx] = max_bin
                # Find which PRS type has this max bin
                max_prs_type = bins_df.iloc[idx].idxmax()
                # Get the corresponding scaled PRS value
                prs_col_name = f'scaled_prs_{max_prs_type}'
                scaled_combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name]
                combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name.replace('scaled_','')]
                n_any_high_risk += 1
                
            elif prop_case_bins.iloc[idx] > 0.6:
                # Priority 2: If >60% case bins (but not any_high_risk), use MAX bin
                max_bin = bins_df.iloc[idx].max()
                combined_bin[idx] = max_bin
                # Find which PRS type has this max bin
                max_prs_type = bins_df.iloc[idx].idxmax()
                # Get the corresponding scaled PRS value
                prs_col_name = f'scaled_prs_{max_prs_type}'
                scaled_combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name]
                combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name.replace('scaled_','')]
                n_prop_high += 1
                
            else:
                # If any_high_risk = 0, use MIN bin
                min_bin = bins_df.iloc[idx].min()
                combined_bin[idx] = min_bin
                # Find which PRS type has this min bin
                min_prs_type = bins_df.iloc[idx].idxmin()
                # Get the corresponding scaled PRS value
                prs_col_name = f'scaled_prs_{min_prs_type}'
                scaled_combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name]
                combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name.replace('scaled_','')]
                n_prop_low += 1
#       for idx in range(len(holdout_processed)):
#           if valid_prs_count.iloc[idx] == 0:
#               combined_bin[idx] = np.nan
#               combined_prs[idx] = np.nan
#               scaled_combined_prs[idx] = np.nan
#           elif prop_case_bins.iloc[idx] > 0.6:
#               # If >60% case bins, use MAX bin
#               max_bin = bins_df.iloc[idx].max()
#               combined_bin[idx] = max_bin
#               # Find which PRS type has this max bin
#               max_prs_type = bins_df.iloc[idx].idxmax()
#               # Get the corresponding scaled PRS value
#               prs_col_name = f'scaled_prs_{max_prs_type}'
#               scaled_combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name]
#               combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name.replace('scaled_','')]
#               
#           else:
#               # If ≤60% case bins, use MIN bin
#               min_bin = bins_df.iloc[idx].min()
#               combined_bin[idx] = min_bin
#               # Find which PRS type has this min bin
#               min_prs_type = bins_df.iloc[idx].idxmin()
#               # Get the corresponding scaled PRS value
#               prs_col_name = f'scaled_prs_{min_prs_type}'
#               scaled_combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name]
#               combined_prs[idx] = holdout_processed.loc[holdout_processed.index[idx], prs_col_name.replace('scaled_','')]
                
                
        holdout_processed['combined_centile_bin'] = combined_bin
        holdout_processed['scaled_prs_combined'] = scaled_combined_prs
        holdout_processed['prs_combined'] = combined_prs
        
        # Assign combined high-risk (bin > threshold_bin)
        holdout_processed['combined_high_risk'] = (combined_bin >= threshold_bin).astype(int)
        
        # Statistics
        combined_high_risk_count = holdout_processed['combined_high_risk'].sum()
        
        print(f"\nAssignment Strategy Breakdown:")
        print(f"{'='*60}")
        print(f"  Priority 1 - any_high_risk > 0 (MAX bin):    {n_any_high_risk:6d} ({n_any_high_risk/len(holdout_processed)*100:5.1f}%)")
        print(f"  Priority 2 - prop_case_bins > 60% (MAX bin): {n_prop_high:6d} ({n_prop_high/len(holdout_processed)*100:5.1f}%)")
        print(f"  Priority 3 - prop_case_bins ≤ 60% (MIN bin): {n_prop_low:6d} ({n_prop_low/len(holdout_processed)*100:5.1f}%)")
        print(f"  Missing data:                                 {n_missing:6d} ({n_missing/len(holdout_processed)*100:5.1f}%)")
        print(f"{'='*60}")
        print(f"Total individuals assigned MAX bin: {n_any_high_risk + n_prop_high} ({(n_any_high_risk + n_prop_high)/len(holdout_processed)*100:.1f}%)")
        print(f"Total individuals assigned MIN bin: {n_prop_low} ({n_prop_low/len(holdout_processed)*100:.1f}%)")
        print(f"\nCombined high-risk (bin >= {threshold_bin}): {combined_high_risk_count}/{len(holdout_processed)} ({combined_high_risk_count/len(holdout_processed)*100:.1f}%)")
        
        # Calculate combined PRS score (mean of PRS values for combined approach)
#       prs_cols_for_combined = [col for col in bin_stats.keys() if col in holdout_processed.columns]
#       if prs_cols_for_combined:
#           holdout_processed['combined_prs'] = holdout_processed[prs_cols_for_combined].mean(axis=1)
            
#       # Assign combined high-risk (bin > threshold_bin)
#       holdout_processed['combined_high_risk'] = (combined_bin >= threshold_bin).astype(int)
#       
#       # Statistics
#       high_case_prop = (prop_case_bins > 0.6).sum()
#       low_case_prop = (prop_case_bins <= 0.6).sum()
#       combined_high_risk_count = holdout_processed['combined_high_risk'].sum()
#       
#       print(f"Individuals with >30% case bins (using MAX): {high_case_prop} ({high_case_prop/len(holdout_processed)*100:.1f}%)")
#       print(f"Individuals with ≤30% case bins (using MIN): {low_case_prop} ({low_case_prop/len(holdout_processed)*100:.1f}%)")
#       print(f"Combined high-risk (bin > {threshold_bin}): {combined_high_risk_count}/{len(holdout_processed)} ({combined_high_risk_count/len(holdout_processed)*100:.1f}%)")
#       print(f"\nCase bin proportion statistics:")
#       print(f"  Mean: {prop_case_bins.mean():.3f}")
#       print(f"  Median: {prop_case_bins.median():.3f}")
#       print(f"  Range: {prop_case_bins.min():.3f} - {prop_case_bins.max():.3f}")
        
    return holdout_processed

def save_bin_statistics(scoresPath,bin_stats, training_stats=None, filename_prefix='training_prs_statistics'):
    """
    Save bin statistics and training statistics to CSV files.
    
    Parameters:
    bin_stats: Dictionary from create_risk_bins()
    training_stats: Dictionary with mean/std for each PRS type
    filename_prefix: Prefix for output files
    """
    
    # Save bin statistics for each PRS type
    for prs_col, stats in bin_stats.items():
        bin_df = stats['bin_dataframe']
        filename = f"{filename_prefix}_{prs_col}_bins.csv"
        bin_df.to_csv(f'{scoresPath}/{filename}', index=False)
        print(f"Bin statistics for {prs_col} saved to {filename}")
        
    # Save training statistics (mean/std) for scaling
    if training_stats is not None:
        training_stats_df = pd.DataFrame.from_dict(training_stats, orient='index')
        training_stats_filename = f"{filename_prefix}_training_stats.csv"
        training_stats_df.to_csv(f'{scoresPath}/{training_stats_filename}')
        print(f"Training statistics saved to {training_stats_filename}")
        
def load_bin_statistics(scoresPath,filename_prefix='training_prs_statistics'):
    """
    Load bin statistics and training statistics from CSV files.
    
    Parameters:
    filename_prefix: Prefix used when saving files
    
    Returns:
    bin_stats: Dictionary of bin statistics
    training_stats: Dictionary of training statistics
    """
    
    prs_columns = ['scaled_prs_main', 'scaled_prs_epi', 'scaled_prs_epi+main', 
                   'scaled_prs_cardio', 'scaled_prs_all']
    
    bin_stats = {}
    
    # Load bin statistics for each PRS type
    for prs_col in prs_columns:
        filename = f"{scoresPath}/{filename_prefix}_{prs_col}_bins.csv"
        try:
            bin_df = pd.read_csv(filename)
            # Reconstruct bin_stats structure
            bin_stats[prs_col] = {
                'bin_dataframe': bin_df,
                'total_bins': len(bin_df)
            }
            print(f"Loaded bin statistics for {prs_col} from {filename}")
        except FileNotFoundError:
            print(f"Warning: {filename} not found, skipping {prs_col}")
            
    # Load training statistics
    training_stats = None
    training_stats_filename = f"{scoresPath}/{filename_prefix}_training_stats.csv"
    try:
        training_stats_df = pd.read_csv(training_stats_filename, index_col=0)
        training_stats = training_stats_df.to_dict('index')
        print(f"Loaded training statistics from {training_stats_filename}")
    except FileNotFoundError:
        print(f"Warning: {training_stats_filename} not found")
        
    return bin_stats, training_stats


# Example usage and main execution
if __name__ == "__main__":
    
    pheno='type2Diabetes'
    scoresPath = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi/scores'
    figPath = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi/figures'
    trainingPath = f'{scoresPath}/combinedPRSGroups.csv'
    holdoutPath = f'{scoresPath}/combinedPRSGroups.holdout.csv'
#   # Example: Create sample data for demonstration
#   print("Creating sample data for demonstration...")
#   np.random.seed(42)
#   
#   # Sample training data
#   n_train = 10000
#   training_data = {
#       'scaled_prs_main': np.random.normal(0, 1, n_train),
#       'scaled_prs_epi': np.random.normal(0, 1, n_train),
#       'scaled_prs_epi+main': np.random.normal(0, 1, n_train),
#       'scaled_prs_cardio': np.random.normal(0, 1, n_train),
#       'scaled_prs_all': np.random.normal(0, 1, n_train)
#   }
#   
#   # Create training DataFrame with IID as index
#   training_df = pd.DataFrame(training_data)
#   training_df.index.name = 'IID'
#   training_df.index = [f'TRAIN_{i}' for i in range(len(training_df))]
#   
#   # Sample holdout data (unscaled PRS)
#   n_holdout = 2000
#   holdout_data = {
#       'prs_main': np.random.normal(2, 3, n_holdout),        # Different mean/std
#       'prs_epi': np.random.normal(-1, 2, n_holdout),
#       'prs_epi+main': np.random.normal(1, 2.5, n_holdout),
#       'prs_cardio': np.random.normal(0.5, 1.5, n_holdout),
#       'prs_all': np.random.normal(-0.5, 2, n_holdout)
#   }
#   
#   # Create holdout DataFrame with IID as index
#   holdout_df = pd.DataFrame(holdout_data)
#   holdout_df.index.name = 'IID'
#   holdout_df.index = [f'HOLDOUT_{i}' for i in range(len(holdout_df))]
    
    #download combinedPRS training Data
    trainingDf = pd.read_csv(trainingPath)
    
    #see if the analysis has been done so as not to overwrite raw file
    try:
        holdoutDf = pd.read_csv(f'{scoresPath}/combinedPRSGroups.holdoutRaw.csv')
        holdoutDf.to_csv(holdoutPath,index=False)
    except FileNotFoundError:
        holdoutDf = pd.read_csv(holdoutPath)
        
        
    try:	
        trainingDf.set_index('IID',inplace=True)
    except ValueError: #the IID is an index so needs to be reindexed
        trainingDf.rename(columns={'Unnamed: 0':'IID'},inplace=True)
        trainingDf.set_index('IID',inplace=True)
        
    try:	
        holdoutDf.rename(columns={'Unnamed: 0':'IID'},inplace=True)
    except KeyError: #the IID is not an index
        pass

    
    # Step 1: Create risk bins from training data
    bin_stats = create_risk_bins(trainingDf, n_bins=1000)
    
    # Step 2: Calculate training statistics for manual scaling
#   training_stats = calculate_training_statistics(trainingDf)
    
    # Step 3: Assign risk thresholds to holdout data
    holdout_processed = assign_risk_thresholds(
        holdoutDf, 
        bin_stats, 
        training_scalers=None,  # Use training stats for manual scaling
        threshold_percentile=80
    )
    
    # Step 4: Save results to CSV files
    save_bin_statistics(scoresPath,bin_stats, filename_prefix='prs_statistics')
    
    holdoutDf.to_csv(f'{scoresPath}/combinedPRSGroups.holdoutRaw.csv',index=False)
    holdout_processed.to_csv(holdoutPath)
    
    #get an odds ratio for each scaled_prs
    prs_cols = [col for col in holdout_processed.columns if 'scaled_prs' in col and 'threshold' not in col]
    percentileORPRS = calculate_odds_ratio_for_prs(holdout_processed[prs_cols+['PHENOTYPE',"IID"]],'scaled_prs')
    percentileORBin = calculate_odds_ratio_for_prs(holdout_processed[['IID','combined_centile_bin','PHENOTYPE']],'centile_bin')
    
    
    #match file names to validation set
    outputORFile = f'{scoresPath}/combinedORPRSGroups.holdout.csv'
    
    percentileOR = pd.concat([percentileORPRS,percentileORBin],ignore_index=True)
    percentileOR.to_csv(outputORFile,index=False)
    
    prsDf = holdout_processed[['IID','combined_centile_bin','PHENOTYPE']].rename(columns={'combined_centile_bin':'scaled_prs'})
    create_prs_plots(prsDf,'combined',figPath,scoresPath,'mixed.holdout')
    
    create_case_control_histogram(holdout_processed,'PHENOTYPE','combined_centile_bin',figPath,figsize=(12,6))
    
    # Display sample results
    print("\n" + "="*60)
    print("SAMPLE RESULTS")
    print("="*60)
    
    print("\nBin statistics for scaled_prs_main (first 5 bins):")
    if 'scaled_prs_main' in bin_stats:
        print(bin_stats['scaled_prs_main']['bin_dataframe'].head())
        
#   print("\nTraining statistics:")
#   if training_stats:
#       for prs_type, stats in training_stats.items():
#           print(f"  {prs_type}: mean={stats['mean']:.4f}, std={stats['std']:.4f}")
            
    print(f"\nHoldout data with risk assignments (first 10 rows):")
    display_cols = [col for col in holdout_processed.columns if 'high_risk' in col or 'threshold' in col][:5]
    print(holdout_processed[display_cols].head(10))
    
    print(f"\nSummary of high-risk assignments:")
    for col in holdout_processed.columns:
        if col.endswith('_high_risk'):
            count = holdout_processed[col].sum()
            pct = count / len(holdout_processed) * 100
            print(f"  {col}: {count}/{len(holdout_processed)} ({pct:.1f}%)")