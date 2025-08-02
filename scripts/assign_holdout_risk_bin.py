#!/usr/bin/env python3

import pandas as pd
import numpy as np

def create_risk_bins(training_df, n_bins=1000):
    """
    Create risk bins for PRS scores and calculate bin statistics.
    
    Parameters:
    training_df: DataFrame with PRS columns and IID as index
    n_bins: Number of risk bins to create (default 1000)
    
    Returns:
    bin_stats: Dictionary containing bin statistics for each PRS type
    """
    
    prs_columns = ['scaled_prs_main', 'scaled_prs_epi', 'scaled_prs_epi+main', 
                   'scaled_prs_cardio', 'scaled_prs_all']
    
    bin_stats = {}
    
    print(f"Creating {n_bins} risk bins for each PRS type...")
    print(f"Training data shape: {training_df.shape}")
    
    for prs_col in prs_columns:
        if prs_col not in training_df.columns:
            print(f"Warning: {prs_col} not found in training data. Skipping...")
            continue
        
        print(f"\nProcessing {prs_col}...")
        
        # Create bins using quantiles to ensure equal sample sizes
        prs_values = training_df[prs_col].dropna()
        
        # Create bin edges using quantiles
        bin_edges = np.linspace(0, 1, n_bins + 1)
        quantile_edges = prs_values.quantile(bin_edges).values
        
        # Assign bin numbers (0 to n_bins-1)
        bin_assignments = pd.cut(prs_values, bins=quantile_edges, 
                               labels=range(n_bins), include_lowest=True)
        
        # Calculate bin statistics
        bin_data = []
        for bin_num in range(n_bins):
            bin_mask = bin_assignments == bin_num
            bin_values = prs_values[bin_mask]
            
            if len(bin_values) > 0:
                bin_stats_dict = {
                    'bin_number': bin_num,
                    'min_value': bin_values.min(),
                    'max_value': bin_values.max(),
                    'mean_value': bin_values.mean(),
                    'count': len(bin_values),
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
                    'count': 0,
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
            'total_bins': n_bins
        }
        
        print(f"  Created {n_bins} bins")
        print(f"  Value range: {prs_values.min():.4f} to {prs_values.max():.4f}")
        print(f"  Non-empty bins: {bin_df['count'].gt(0).sum()}")
        
    return bin_stats

def assign_risk_thresholds(holdout_df, bin_stats, training_scalers=None, threshold_percentile=80):
    """
    Assign risk thresholds to holdout data based on training bin statistics.
    
    Parameters:
    holdout_df: DataFrame with PRS columns (may be unscaled)
    bin_stats: Dictionary from create_risk_bins()
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
    
    # Assign bins and thresholds for each PRS type
    for prs_col, stats in bin_stats.items():
        if prs_col not in holdout_processed.columns:
            print(f"Warning: {prs_col} not found in holdout data. Skipping...")
            continue
        
        print(f"\nProcessing {prs_col} in holdout data...")
        
        bin_df = stats['bin_dataframe']
        bin_edges = stats['bin_edges']
        
        # Assign bins to holdout data
        prs_values = holdout_processed[prs_col]
        bin_assignments = pd.cut(prs_values, bins=bin_edges, 
                               labels=range(len(bin_edges)-1), include_lowest=True)
        
        # Convert to numeric (handling NaN)
        bin_assignments_numeric = pd.to_numeric(bin_assignments, errors='coerce')
        
        # Create column names for this PRS type
        bin_col = f'{prs_col}_bin'
        threshold_col = f'{prs_col}_high_risk'
        
        # Assign bin numbers
        holdout_processed[bin_col] = bin_assignments_numeric
        
        # Assign high-risk status (1 if bin >= threshold_bin, 0 otherwise)
        holdout_processed[threshold_col] = (bin_assignments_numeric >= threshold_bin).astype(int)
        
        # Add threshold value for reference
        threshold_value = bin_df.loc[threshold_bin, 'min_value'] if threshold_bin < len(bin_df) else bin_df['max_value'].max()
        holdout_processed[f'{prs_col}_threshold_value'] = threshold_value
        
        # Statistics
        high_risk_count = holdout_processed[threshold_col].sum()
        total_count = len(holdout_processed.dropna(subset=[prs_col]))
        
        print(f"  Assigned bins: 0-{bin_assignments_numeric.max()}")
        print(f"  High-risk threshold value: {threshold_value:.4f}")
        print(f"  High-risk individuals: {high_risk_count}/{total_count} ({high_risk_count/total_count*100:.1f}%)")
        
    # Create composite high-risk indicator (high-risk in ANY PRS type)
    risk_columns = [col for col in holdout_processed.columns if col.endswith('_high_risk')]
    if risk_columns:
        holdout_processed['any_high_risk'] = holdout_processed[risk_columns].max(axis=1)
        any_high_risk_count = holdout_processed['any_high_risk'].sum()
        print(f"\nComposite high-risk (any PRS type): {any_high_risk_count}/{len(holdout_processed)} ({any_high_risk_count/len(holdout_processed)*100:.1f}%)")
        
    return holdout_processed

def save_bin_statistics(bin_stats, training_stats=None, filename_prefix='prs_statistics'):
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
        bin_df.to_csv(filename, index=False)
        print(f"Bin statistics for {prs_col} saved to {filename}")
        
    # Save training statistics (mean/std) for scaling
    if training_stats is not None:
        training_stats_df = pd.DataFrame.from_dict(training_stats, orient='index')
        training_stats_filename = f"{filename_prefix}_training_stats.csv"
        training_stats_df.to_csv(training_stats_filename)
        print(f"Training statistics saved to {training_stats_filename}")
        
def load_bin_statistics(filename_prefix='prs_statistics'):
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
        filename = f"{filename_prefix}_{prs_col}_bins.csv"
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
    training_stats_filename = f"{filename_prefix}_training_stats.csv"
    try:
        training_stats_df = pd.read_csv(training_stats_filename, index_col=0)
        training_stats = training_stats_df.to_dict('index')
        print(f"Loaded training statistics from {training_stats_filename}")
    except FileNotFoundError:
        print(f"Warning: {training_stats_filename} not found")
        
    return bin_stats, training_stats

def calculate_training_statistics(training_df):
    """
    Calculate mean and standard deviation for each PRS column in training data.
    
    Parameters:
    training_df: Training DataFrame with scaled PRS columns
    
    Returns:
    training_stats: Dictionary with mean and std for each PRS type
    """
    
    prs_columns = ['scaled_prs_main', 'scaled_prs_epi', 'scaled_prs_epi+main', 
                   'scaled_prs_cardio', 'scaled_prs_all']
    
    training_stats = {}
    
    print("Calculating training statistics for scaling...")
    
    for prs_col in prs_columns:
        if prs_col in training_df.columns:
            prs_values = training_df[prs_col].dropna()
            
            stats = {
                'mean': float(prs_values.mean()),
                'std': float(prs_values.std()),
                'count': len(prs_values),
                'min': float(prs_values.min()),
                'max': float(prs_values.max())
            }
            
            training_stats[prs_col] = stats
            
            print(f"  {prs_col}:")
            print(f"    Mean: {stats['mean']:.6f}")
            print(f"    Std: {stats['std']:.6f}")
            print(f"    Count: {stats['count']}")
        else:
            print(f"Warning: {prs_col} not found in training data")
            
    return training_stats

def scale_holdout_data_manually(holdout_df, training_stats):
    """
    Manually scale holdout data using training mean and standard deviation.
    
    Parameters:
    holdout_df: DataFrame with unscaled PRS columns
    training_stats: Dictionary with mean/std from calculate_training_statistics()
    
    Returns:
    holdout_scaled: DataFrame with scaled PRS columns added
    """
    
    # Base PRS column names (what we expect in holdout data)
    base_prs_columns = ['prs_main', 'prs_epi', 'prs_epi+main', 'prs_cardio', 'prs_all']
    
    holdout_scaled = holdout_df.copy()
    
    print("Manually scaling holdout data using training statistics...")
    
    for base_col in base_prs_columns:
        scaled_col = f'scaled_{base_col}'
        
        if base_col in holdout_df.columns and scaled_col in training_stats:
            
            # Get training statistics
            train_mean = training_stats[scaled_col]['mean']
            train_std = training_stats[scaled_col]['std']
            
            # Manual scaling: (x - mean) / std
            holdout_scaled[scaled_col] = (holdout_df[base_col] - train_mean) / train_std
            
            # Calculate statistics for verification
            scaled_values = holdout_scaled[scaled_col].dropna()
            
            print(f"  {base_col} -> {scaled_col}")
            print(f"    Training mean: {train_mean:.6f}, std: {train_std:.6f}")
            print(f"    Holdout scaled mean: {scaled_values.mean():.6f}, std: {scaled_values.std():.6f}")
            print(f"    Holdout scaled range: [{scaled_values.min():.3f}, {scaled_values.max():.3f}]")
            
        elif base_col not in holdout_df.columns:
            print(f"Warning: {base_col} not found in holdout data")
        elif scaled_col not in training_stats:
            print(f"Warning: No training statistics found for {scaled_col}")
            
    return holdout_scaled

# Example usage and main execution
if __name__ == "__main__":
    
    # Example: Create sample data for demonstration
    print("Creating sample data for demonstration...")
    np.random.seed(42)
    
    # Sample training data
    n_train = 10000
    training_data = {
        'scaled_prs_main': np.random.normal(0, 1, n_train),
        'scaled_prs_epi': np.random.normal(0, 1, n_train),
        'scaled_prs_epi+main': np.random.normal(0, 1, n_train),
        'scaled_prs_cardio': np.random.normal(0, 1, n_train),
        'scaled_prs_all': np.random.normal(0, 1, n_train)
    }
    
    # Create training DataFrame with IID as index
    training_df = pd.DataFrame(training_data)
    training_df.index.name = 'IID'
    training_df.index = [f'TRAIN_{i}' for i in range(len(training_df))]
    
    # Sample holdout data (unscaled PRS)
    n_holdout = 2000
    holdout_data = {
        'prs_main': np.random.normal(2, 3, n_holdout),        # Different mean/std
        'prs_epi': np.random.normal(-1, 2, n_holdout),
        'prs_epi+main': np.random.normal(1, 2.5, n_holdout),
        'prs_cardio': np.random.normal(0.5, 1.5, n_holdout),
        'prs_all': np.random.normal(-0.5, 2, n_holdout)
    }
    
    # Create holdout DataFrame with IID as index
    holdout_df = pd.DataFrame(holdout_data)
    holdout_df.index.name = 'IID'
    holdout_df.index = [f'HOLDOUT_{i}' for i in range(len(holdout_df))]
    
    # Step 1: Create risk bins from training data
    bin_stats = create_risk_bins(training_df, n_bins=1000)
    
    # Step 2: Calculate training statistics for manual scaling
    training_stats = calculate_training_statistics(training_df)
    
    # Step 3: Assign risk thresholds to holdout data
    holdout_processed = assign_risk_thresholds(
        holdout_df, 
        bin_stats, 
        training_stats=training_stats,  # Use training stats for manual scaling
        threshold_percentile=80
    )
    
    # Step 4: Save results to CSV files
    save_bin_statistics(bin_stats, training_stats, filename_prefix='prs_statistics')
    holdout_processed.to_csv('holdout_with_risk_assignments.csv')
    
    # Display sample results
    print("\n" + "="*60)
    print("SAMPLE RESULTS")
    print("="*60)
    
    print("\nBin statistics for scaled_prs_main (first 5 bins):")
    if 'scaled_prs_main' in bin_stats:
        print(bin_stats['scaled_prs_main']['bin_dataframe'].head())
        
    print("\nTraining statistics:")
    if training_stats:
        for prs_type, stats in training_stats.items():
            print(f"  {prs_type}: mean={stats['mean']:.4f}, std={stats['std']:.4f}")
            
    print(f"\nHoldout data with risk assignments (first 10 rows):")
    display_cols = [col for col in holdout_processed.columns if 'high_risk' in col or 'threshold' in col][:5]
    print(holdout_processed[display_cols].head(10))
    
    print(f"\nSummary of high-risk assignments:")
    for col in holdout_processed.columns:
        if col.endswith('_high_risk'):
            count = holdout_processed[col].sum()
            pct = count / len(holdout_processed) * 100
            print(f"  {col}: {count}/{len(holdout_processed)} ({pct:.1f}%)")