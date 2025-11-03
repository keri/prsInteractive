#!/usr/bin/env python3

#Keri Multerer April 7, 2025
#methods to calculate PRS performance against clinical measures

import numpy as np
import json
import pandas as pd
from sklearn.metrics import roc_curve, auc, roc_auc_score, brier_score_loss
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import calibration_curve
from scipy import stats
import statsmodels.api as sm
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression



def determine_risk_direction(df, measure, outcome_column=None, correlation_threshold=0.1):
    """
    Determine whether higher or lower values of a clinical measure indicate higher risk.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing clinical measures
    measure : str
        Column name of the clinical measure to analyze
    outcome_column : str or None
        If provided, the binary outcome column to correlate with (1=disease/event)
        If None, will use statistical properties of the distribution
    correlation_threshold : float
        Threshold for correlation strength to determine direction
        
    Returns:
    --------
    bool
        True if higher values indicate higher risk, False if lower values indicate higher risk
    """
    if outcome_column is not None and outcome_column in df.columns:
        # If we have an outcome column, use correlation
        if set(df[outcome_column].unique()).issubset({0, 1}):
            # Use point-biserial correlation (equivalent to Pearson with binary variable)
            correlation = df[measure].corr(df[outcome_column])
            return correlation > 0  # Positive correlation means higher values = higher risk
        
    # If no outcome column or correlation is weak, use distributional properties
    # This is a heuristic approach that can be refined based on domain knowledge
        
    # For most clinical measures, positive skew often indicates that higher values are abnormal
    skewness = stats.skew(df[measure].dropna())
    
    # Check if the distribution is notably skewed
    if abs(skewness) > 0.5:
        return skewness > 0  # If positively skewed, higher values likely indicate risk
    
    # If skewness is inconclusive, default to higher values = higher risk
    # This is a reasonable default for many clinical measures
    return True


def convert_to_binary(df, df_validation, clinical_measures, thresholds=None, high_risk_quintile=True, 
                      risk_directions=None, outcome_column=None):
    """
    Convert continuous clinical measures to binary (0/1) based on specified thresholds
    or high/low risk quintile.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing clinical measures
    clinical_measures : list or dict
        If list: column names of clinical measures to convert
        If dict: {measure_name: threshold_value} to use specific thresholds
    thresholds : dict or None
        Dictionary mapping clinical measure names to threshold values
    high_risk_quintile : bool or float
        If True, use quintile (80th or 20th percentile) based on risk direction
        If float between 0 and 1, use that percentile instead of 0.8/0.2
    risk_directions : dict or None
        Dictionary mapping measure names to booleans:
        True if higher values indicate higher risk, False if lower values indicate higher risk
        If None, direction will be determined automatically
    outcome_column : str or None
        Binary outcome column to use for determining risk direction
        
    Returns:
    --------
    pandas.DataFrame
        Copy of input DataFrame with additional binary columns named '{measure}_binary'
    dict
        Dictionary of thresholds used for each measure
    dict
        Dictionary of risk directions determined for each measure
    """
    df_result = df.copy()
    df_result_validation = df_validation.copy()
    used_thresholds = {}
    determined_directions = {}
    
    # Set quintile threshold (default to 0.8 or 0.2 if True)
    if high_risk_quintile is True:
        quintile_threshold = 0.2
    elif isinstance(high_risk_quintile, (int, float)) and 0 <= high_risk_quintile <= 1:
        quintile_threshold = high_risk_quintile
    else:
        quintile_threshold = 0.2  # Default
        
    # Convert clinical_measures list to dict with None values if it's a list
    if isinstance(clinical_measures, list):
        clinical_measures_dict = {measure: None for measure in clinical_measures}
    else:
        clinical_measures_dict = clinical_measures
        
    # Override with provided thresholds if any
    if thresholds:
        for measure, threshold in thresholds.items():
            if measure in clinical_measures_dict:
                clinical_measures_dict[measure] = threshold
                
    # Initialize risk_directions if not provided
    if risk_directions is None:
        risk_directions = {}
        
    for measure, threshold in clinical_measures_dict.items():
        # Skip if the measure doesn't exist in the dataframe
        if measure not in df_result.columns:
            print(f"Warning: Measure '{measure}' not found in dataframe")
            continue
        
        # Check if the measure is already binary
        unique_values = set(df_result[measure].unique())
        if unique_values.issubset({0, 1}) and len(unique_values) <= 2:
            print(f"Note: Measure '{measure}' is already binary, keeping as is")
            df_result[f"{measure}_binary"] = df_result[measure]
            used_thresholds[measure] = "already binary"
            determined_directions[measure] = True  # Default for binary
            continue
        
        # Determine risk direction if not specified
        if measure not in risk_directions:
            higher_is_riskier = determine_risk_direction(df_result, measure, outcome_column)
            determined_directions[measure] = higher_is_riskier
        else:
            higher_is_riskier = risk_directions[measure]
            determined_directions[measure] = higher_is_riskier
            
        # If no threshold is provided and high_risk_quintile is True, use quintile
        if threshold is None and high_risk_quintile:
            # Choose appropriate percentile based on risk direction
            if higher_is_riskier:
                threshold = df_result[measure].quantile(quintile_threshold)  # Default 20th percentile
                direction_str = "higher values → higher risk"
            else:
                threshold = df_result[measure].quantile(1 - quintile_threshold)  # Default 80th percentile
                direction_str = "lower values → higher risk"
                
            print(f"Using {direction_str} for '{measure}', threshold = {threshold:.3f}")
            
        # Apply threshold to create binary variable
        if threshold is not None:
            if higher_is_riskier:
                df_result[f"{measure}_binary"] = (df_result[measure] <= threshold).astype(int)
                df_result_validation[f"{measure}_binary"] = (df_result_validation[measure] <= threshold).astype(int)
            else:
                df_result[f"{measure}_binary"] = (df_result[measure] >= threshold).astype(int)
                df_result_validation[f"{measure}_binary"] = (df_result_validation[measure] >= threshold).astype(int)
                
            used_thresholds[measure] = threshold
        else:
            print(f"Warning: No threshold specified for '{measure}' and quintile option not used")
            
    return df_result, df_result_validation, used_thresholds, determined_directions


def calculate_auc(dfFull, prs_methods, clinical_measures):
    """
    Calculate Area Under the ROC Curve (AUC) for different PRS methods against clinical measures.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_methods : list
        List of column names for the different PRS methods
    clinical_measures : list
        List of column names for the different clinical outcome measures
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing AUC scores for each PRS method against each clinical measure
    """
    resultsAllPRS = pd.DataFrame(index=prs_methods, columns=clinical_measures)
    resultsTop20PRS = pd.DataFrame(index=prs_methods, columns=clinical_measures)
    

    for prs in prs_methods:
        for measure in clinical_measures:
            #get the low risk clinical measure
            df = dfFull[dfFull[measure] == 1]
            
            #get the high prs 
            threshold = df[prs].quantile(.80)
            df['high_prs'] = (df[prs] <= threshold).astype(int)
            
            # Ensure the clinical measure is binary (0/1)
            if not set(df['PHENOTYPE'].unique()).issubset({0, 1}):
                raise ValueError(f"PHENOTYPE must be binary (0/1)")
                
            # Calculate AUC
            try:
                auc = roc_auc_score(df['PHENOTYPE'], df[prs])
                resultsAllPRS.loc[prs, measure] = auc
                
                dfHigh = df[df['high_prs'] == 1]
                auc = roc_auc_score(dfHigh['PHENOTYPE'], dfHigh[prs])
                resultsTop20PRS.loc[prs, measure] = auc
                
            except Exception as e:
                resultsAllPRS.loc[prs, measure] = np.nan
                print(f"Error calculating AUC for {prs} vs low risk {measure}: {e}")
                
                resultsTop20PRS.loc[prs, measure] = np.nan
                print(f"Error calculating AUC for high risk {prs} vs low risk {measure}: {e}")
                
    return resultsAllPRS, resultsTop20PRS

def calculate_single_prs_nri(df,prs_method1,clinical_measure):
    """
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_method1 : str
        Column name for the first PRS method (reference)
    clinical_measure : str
        Column name for the clinical outcome measure (binary, 0/1)

    
    Returns:
    --------
    tuple
        (NRI, NRI for events, NRI for non-events)
    """
    if clinical_measure == 'None':
        df_clinic = df.copy()
        
    else:
        # Filter to low clinical risk only
        df_clinic = df[df[clinical_measure] == 1].copy()
        
    # Categorize PRS
    risk_threshold = df_clinic[prs_method1].quantile(0.8)
    df_clinic['high_prs'] = df_clinic[prs_method1].apply(lambda x: 1 if x >= risk_threshold else 0)
    
    # Split into cases and controls
    events = df_clinic[df_clinic['PHENOTYPE'] == 1]
    non_events = df_clinic[df_clinic['PHENOTYPE'] == 0]
    
    # NRI = (percent of cases reclassified up) - (percent of controls reclassified up)
    nri_cases = events["high_prs"].mean()
    nri_controls = non_events["high_prs"].mean()
    
    prs_nri = nri_cases - nri_controls
    
    return prs_nri, nri_cases, nri_controls

def calculate_nri(df, prs_method1, prs_method2, clinical_measure,clinical):
    """
    Calculate Net Reclassification Improvement (NRI) comparing two PRS methods.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_method1 : str
        Column name for the first PRS method (reference)
    prs_method2 : str
        Column name for the second PRS method (new)
    clinical_measure : str
        Column name for the clinical outcome measure (binary, 0/1)
    clinical : bool
        True is clinical measure is compared to PRS
        
    Returns:
    --------
    tuple
        (NRI, NRI for events, NRI for non-events)
    """
    
    if clinical_measure != "None":
        # Ensure the clinical measure is binary
        if not set(df[clinical_measure].unique()).issubset({0, 1}):
            raise ValueError(f"Clinical measure {clinical_measure} must be binary (0/1)")
            
        # Extract events and non-events
        events = df[clinical_measure] == 1
        non_events = df[clinical_measure] == 0
        
    else:
        
        # Extract events and non-events
        events = df['PHENOTYPE'] == 1
        non_events = df['PHENOTYPE'] == 0
        
    # Create risk categories based on thresholds
    def categorize(scores, thresholds):
        categories = np.zeros(len(scores), dtype=int)
        for i, threshold in enumerate(sorted(thresholds)):
            categories[scores >= threshold] = i + 1
        return categories
    
    # Calculate NRI for all pairs of PRS methods
    if clinical:           
        #find the clinical measure with _binary in label
        if '_binary' in prs_method1:
            # Default to quartiles of the first PRS method
            risk_threshold = df[prs_method2].quantile(0.8)
            # Categorize: Top 20% as "High Risk", others as "Low Risk"
            risk_cat2 = df[prs_method2].apply(lambda x: 1 if x >= risk_threshold else 0)

            risk_cat1 = df[prs_method1].apply(lambda x: 1 if x == 0 else 0)
        
        else:
#           risk_thresholds = list(df[prs_method1].quantile([0.20, 0.5, 0.80]))
#           risk_cat2 = categorize(df[prs_method1], risk_thresholds)
#           risk_cat1 = df[prs_method2].astype('category')
            risk_threshold = df[prs_method1].quantile(0.8)
            # Categorize: Top 20% as "High Risk", others as "Low Risk"
            risk_cat1 = df[prs_method1].apply(lambda x: 1 if x >= risk_threshold else 0)
            
            risk_cat2 = df[prs_method2].apply(lambda x: 1 if x == 0 else 0)
            
    else:
        # Default to quartiles of the first PRS method
#       risk_thresholds = list(df[prs_method2].quantile([0.20, 0.5, 0.8]))
#       risk_cat1 = categorize(df[prs_method2], risk_thresholds)
#       
#       risk_thresholds = list(df[prs_method1].quantile([0.20, 0.5, 0.80]))
#       risk_cat2 = categorize(df[prs_method1], risk_thresholds)
        risk_threshold = df[prs_method1].quantile(0.8)
        # Categorize: Top 20% as "High Risk", others as "Low Risk"
        risk_cat1 = df[prs_method1].apply(lambda x: 1 if x >= risk_threshold else 0)
        
        risk_threshold = df[prs_method2].quantile(0.8)
        # Categorize: Top 20% as "High Risk", others as "Low Risk"
        risk_cat2 = df[prs_method2].apply(lambda x: 1 if x >= risk_threshold else 0)
        

    
    # Calculate proportions of up-classifications and down-classifications
    up_events = np.mean(risk_cat2[events] > risk_cat1[events])
    down_events = np.mean(risk_cat2[events] < risk_cat1[events])
    
    up_non_events = np.mean(risk_cat2[non_events] > risk_cat1[non_events])
    down_non_events = np.mean(risk_cat2[non_events] < risk_cat1[non_events])
    
    # Calculate NRI components
    nri_events = up_events - down_events
    nri_non_events = down_non_events - up_non_events
    nri = nri_events + nri_non_events
        
    return nri, nri_events, nri_non_events

def calculate_calibration_comparing_prs(df,prs_method1,prs_method2,n_bins=10):
    """
    Calculate calibration metrics for different PRS methods compared to eachother.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_method1 and prs_method2: strings
        prs_methods to compare
    n_bins : int
        Number of bins for calibration curve
        
    Returns:
    --------
    dict
        Dictionary containing calibration metrics for each PRS method against each clinical measure
    """
    
    #create a binary value for people that are in high risk for prs1
    df["high_prs1"] = (df[prs_method1] >= .8).astype(int)
    
    try:
        
           # Normalize PRS to 0-1 scale if needed
        if df[prs_method2].min() < 0 or df[prs_method2].max() > 1:
            # Min-Max scaling
            prs_scaled = (df[prs_method2] - df[prs_method2].min()) / (df[prs_method2].max() - df[prs_method2].min())
        else:
            prs_scaled = df[prs_method2]
            
        # Calculate calibration curve
        prob_true, prob_pred = calibration_curve(df["high_prs1"], prs_scaled, n_bins=n_bins)
        
        # Calculate Brier score
        brier = brier_score_loss(df["high_prs1"], prs_scaled)
        
        # Calculate Hosmer-Lemeshow statistic
        # Group predictions into bins
        quantiles = np.percentile(prs_scaled, np.linspace(0, 100, n_bins+1))
        bins = np.digitize(prs_scaled, quantiles)
        
        observed = np.zeros(n_bins)
        expected = np.zeros(n_bins)
        counts = np.zeros(n_bins)
        
        for i in range(1, n_bins+1):
            bin_indices = (bins == i)
            if np.sum(bin_indices) > 0:
                observed[i-1] = np.sum(df["high_prs1"][bin_indices])
                expected[i-1] = np.sum(prs_scaled[bin_indices])
                counts[i-1] = np.sum(bin_indices)
                
        # Calculate chi-square statistic
        with np.errstate(divide='ignore', invalid='ignore'):
            chi_square = np.nansum((observed - expected)**2 / (expected * (1 - expected/counts)))
            
        p_value = 1 - stats.chi2.cdf(chi_square, n_bins-2)
        
        # Store results
        results = {
            'prob_true': prob_true,
            'prob_pred': prob_pred,
            'brier_score': brier,
            'hosmer_lemeshow_chi2': chi_square,
            'hosmer_lemeshow_pvalue': p_value
        }
        
    except Exception as e:
        results = {
            'error': str(e)
        }
        print(f"Error in calibration for {prs_method1} vs {prs_method2}: {e}")
        
    return results
    
    


def calculate_calibration(dfFull, prs_methods, clinical_measures, pos_label, n_bins=10):
    """
    Calculate calibration metrics for different PRS methods against clinical measures.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_methods : list
        List of column names for the different PRS methods
    clinical_measures : list
        List of column names for the different clinical outcome measures
    n_bins : int
        Number of bins for calibration curve
        
    Returns:
    --------
    dict
        Dictionary containing calibration metrics for each PRS method against each clinical measure
    """
    results = {}
    auc_plots = {}
    for prs in prs_methods:
        results[prs] = {}
        auc_plots[prs] = {}

        for measure in clinical_measures+['None']:
            
            if measure != 'None':
                
                #get the people with low clinical measures
                df = dfFull[dfFull[measure] == 1]
    
                # Ensure the clinical measure is binary
                if not set(df[measure].unique()).issubset({0, 1}):
                    print(f"Skipping {measure} as it's not binary")
                    continue
            
            else:
                df = dfFull.copy()
                
            #recalibrate PRS for people with low clincal measures
            
            # Use your PRS score as "prediction" (raw, uncalibrated) and actual phenotype
            y_true = df['PHENOTYPE']   # 0/1 outcome
            prs_preds = df[prs] # change to your PRS column
            
#           # Step 2: Fit Isotonic Regression
#           iso = IsotonicRegression(out_of_bounds='clip')
#           prs_calibrated = iso.fit_transform(prs_preds, y_true)
#           
#           # Save calibrated predictions
#           df['prs_calibrated'] = prs_calibrated
#           
            #           # Normalize PRS to 0-1 scale if needed
            if df[prs].min() < 0 or df[prs].max() > 1:
                # Min-Max scaling
                prs_scaled = (df[prs] - df[prs].min()) / (df[prs].max() - df[prs].min())
            else:
                prs_scaled = df[prs]
                            
            try:
                # Calculate calibration curve
                prob_true, prob_pred = calibration_curve(df[pos_label], prs_scaled, n_bins=n_bins)
                
                # Calculate Brier score
                brier = brier_score_loss(df[pos_label], prs_scaled)
                
                # Calculate Hosmer-Lemeshow statistic
                # Group predictions into bins
                quantiles = np.percentile(prs_scaled, np.linspace(0, 100, n_bins+1))
                bins = np.digitize(prs_scaled, quantiles)
                
                observed = np.zeros(n_bins)
                expected = np.zeros(n_bins)
                counts = np.zeros(n_bins)
                
                for i in range(1, n_bins+1):
                    bin_indices = (bins == i)
                    if np.sum(bin_indices) > 0:
                        observed[i-1] = np.sum(df[pos_label][bin_indices])
                        expected[i-1] = np.sum(prs_scaled[bin_indices])
                        counts[i-1] = np.sum(bin_indices)
                        
                # Calculate chi-square statistic
                with np.errstate(divide='ignore', invalid='ignore'):
                    chi_square = np.nansum((observed - expected)**2 / (expected * (1 - expected/counts)))
                    
                p_value = 1 - stats.chi2.cdf(chi_square, n_bins-2)
                
                # Store results
                results[prs][measure] = {
                    'prob_true': prob_true,
                    'prob_pred': prob_pred,
                    'brier_score': brier,
                    'hosmer_lemeshow_chi2': chi_square,
                    'hosmer_lemeshow_pvalue': p_value
                }
                
                #plot auc plots with calibrated PRS
                #auc_plots[prs][measure] = plot_roc_curves(df,'prs_calibrated',measure)
                
            except Exception as e:
                results[prs][measure] = {
                    'error': str(e)
                }
                print(f"Error in calibration for {prs} vs {measure}: {e}")
                
            
        
    return results


#def plot_roc_curves(df, prs_methods, clinical_measure, figsize=(10, 8)):
def plot_roc_curves(df, prs, clinical_measure, figsize=(10, 8)):
    """
    Plot ROC curves comparing different PRS methods for a specific clinical measure.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures
    prs_methods : list
        List of column names for the different PRS methods
    clinical_measure : str
        Column name for the clinical outcome measure
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        ROC curve plot
    """
    plt.figure(figsize=figsize)
    if clinical_measure == 'None':
        clinical_measure = 'PHENOTYPE'
    
#   for prs in prs_methods:
    try:
        fpr, tpr, _ = roc_curve(df[clinical_measure], df[prs])
        auc = roc_auc_score(df[clinical_measure], df[prs])
        plt.plot(fpr, tpr, label=f'{prs} (AUC = {auc:.3f})')
    except Exception as e:
        print(f"Error plotting ROC for {prs}: {e}")
            
    plt.plot([0, 1], [0, 1], 'k--', label='random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curves for {clinical_measure}')
    plt.legend()
    plt.grid(alpha=0.3)
    
    return plt.gcf()


def plot_calibration_curves(calibration_results, prs_methods, clinical_measure, figsize=(10, 8)):
    """
    Plot calibration curves comparing different PRS methods for a specific clinical measure.
    
    Parameters:
    -----------
    calibration_results : dict
        Dictionary from calculate_calibration function
    prs_methods : list
        List of PRS methods to include in the plot
    clinical_measure : str
        Name of the clinical measure
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        Calibration curve plot
    """
    plt.figure(figsize=figsize)
    
    # Plot reference line (perfect calibration)
    plt.plot([0, 1], [0, 1], 'k--', label='Perfectly calibrated')
    
    for prs in prs_methods:
        if clinical_measure in calibration_results[prs] and 'error' not in calibration_results[prs][clinical_measure]:
            result = calibration_results[prs][clinical_measure]
            plt.plot(
                result['prob_pred'], 
                result['prob_true'], 
                'o-', 
                label=f'{prs} (Brier: {result["brier_score"]:.3f})'
            )
            
    plt.xlabel('Predicted probability')
    plt.ylabel('True probability')
    plt.title(f'Calibration Curves for {clinical_measure}')
    plt.legend()
    plt.grid(alpha=0.3)
    
    return plt.gcf()


def plot_risk_distribution(df, measure, binary_col, threshold, higher_is_riskier=True, figsize=(12, 6)):
    """
    Plot the distribution of a clinical measure, highlighting the high-risk threshold.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing clinical measures
    measure : str
        Column name of the continuous clinical measure
    binary_col : str
        Column name of the binarized version of the measure
    threshold : float
        Threshold value used for binarization
    higher_is_riskier : bool
        Whether higher values indicate higher risk
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        Distribution plot
    """
    plt.figure(figsize=figsize)
    
    # Calculate percentile of threshold
    percentile = stats.percentileofscore(df[measure].dropna(), threshold)
    
#   # Plot distribution
#   plt.hist(df[measure], bins=30, alpha=0.6, color='skyblue', edgecolor='black')
    
    # Split data
    cases = df[df['PHENOTYPE'] == 1][measure].dropna()
    controls = df[df['PHENOTYPE'] == 0][measure].dropna()
    
    # Plot histograms
    plt.hist(controls, bins=30, alpha=0.6, label='Controls', color='skyblue', edgecolor='black')
    plt.hist(cases, bins=30, alpha=0.6, label='Cases', color='salmon', edgecolor='black')
    
    # Plot threshold line
    plt.axvline(x=threshold, color='red', linestyle='--', linewidth=2)
    
    # Shaded risk region
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    
    if higher_is_riskier:
        plt.fill_betweenx([0, y_max], threshold, x_max, color='red', alpha=0.1)
        risk_text = f"High Risk\n({100 - percentile:.1f}%)"
        text_x = threshold + (x_max - threshold) * 0.5
    else:
        plt.fill_betweenx([0, y_max], x_min, threshold, color='red', alpha=0.1)        
        risk_text = f"High Risk\n({percentile:.1f}%)"
        text_x = x_min + (threshold - x_min) * 0.5
        
    # Add annotation
    plt.text(text_x, y_max * 0.8, risk_text, ha='center', color='darkred', fontsize=12)
        
    
#   # Plot threshold line
#   plt.axvline(x=threshold, color='red', linestyle='--', linewidth=2)
#   
#   # Determine shaded area (high risk)
#   x_range = plt.xlim()
#   if higher_is_riskier:
#       plt.fill_between([threshold, x_range[1]], 0, plt.ylim()[1], 
#                        color='red', alpha=0.1)
#       risk_text = f"High Risk\n({100-percentile:.1f}%)"
#       text_pos = (threshold + (x_range[1] - threshold) * 0.5, plt.ylim()[1] * 0.8)
#   else:
#       plt.fill_between([x_range[0], threshold], 0, plt.ylim()[1], 
#                        color='red', alpha=0.1)
#       risk_text = f"High Risk\n({percentile:.1f}%)"
#       text_pos = (x_range[0] + (threshold - x_range[0]) * 0.5, plt.ylim()[1] * 0.8)
        
#   plt.text(text_pos[0], text_pos[1], risk_text, 
#            horizontalalignment='center', color='darkred', fontsize=12)
    
    # Calculate percentage of high risk
    high_risk_percent = (df[binary_col] == 1).mean() * 100
    
    plt.title(f'Distribution of {measure}\nThreshold: {threshold:.2f} ({percentile:.1f}th percentile, {high_risk_percent:.1f}% classified as high risk)')
    plt.xlabel(measure)
    plt.ylabel('Frequency')
    plt.grid(alpha=0.3)
    
    return plt.gcf()


#def compare_prs_performance(df, df_validation, prs_methods, clinical_measures, binary_to_use=None, clinical_thresholds=None, risk_thresholds=None, high_risk_quintile=True, risk_directions=None,
#                           outcome_column=None):
#   """
#   Comprehensive comparison of multiple PRS methods against multiple clinical measures
#   
#   Parameters:
#   -----------
#   df : pandas.DataFrame
#       DataFrame containing PRS scores and clinical outcome measures for training data
#   df_validation : pandas.DataFrame 
#       DataFrame containing PRS scores and clinical outcome measures for validation set
#   prs_methods : list
#       List of column names for the different PRS methods
#   clinical_measures : list
#       List of column names for the different clinical outcome measures
#   clinical_thresholds : dict or None
#       Dictionary mapping clinical measure names to threshold values for binarization
#   risk_thresholds : list or None
#       Risk thresholds for NRI calculation. If None, quartiles of the first PRS method will be used.
#   high_risk_quintile : bool or float
#       If True, use quintile (80th or 20th percentile) based on risk direction
#       If float between 0 and 1, use that percentile instead of 0.8/0.2
#   risk_directions : dict or None
#       Dictionary mapping measure names to booleans:
#       True if higher values indicate higher risk, False if lower values indicate higher risk
#   binary_to_use : str or None
#       if string use the validation data to run analysis
#   outcome_column : str or None
#       Binary outcome column to use for determining risk direction
#       
#   Returns:
#   --------
#   dict
#       Dictionary containing performance metrics
#   pandas.DataFrame
#       DataFrame with added binary columns
#   """
#   results = {}
#   
#   # First convert clinical measures to binary if needed
#   df_binary_test, df_binary_validation, used_thresholds, determined_directions = convert_to_binary(
#       df, df_validation, clinical_measures, thresholds=clinical_thresholds, 
#       high_risk_quintile=high_risk_quintile, risk_directions=risk_directions,
#       outcome_column=outcome_column
#   )
#   
#   # Store the thresholds and directions used
#   results['thresholds_used'] = used_thresholds
#   results['risk_directions'] = determined_directions
#   
#   # Create list of binary clinical measures
#   binary_measures = [f"{measure}_binary" for measure in clinical_measures 
#                      if f"{measure}_binary" in df_binary_test.columns]
#   
#   # If no binary measures were created, use original measures (assuming they're already binary)
#   if not binary_measures:
#       binary_measures = clinical_measures
#       
#   if binary_to_use:
#       df_binary = df_binary_validation.copy()
#   else:
#       df_binary = df_binary_test.copy()
#       
#   # Calculate AUC for all methods and binary measures
#   resultsLowClinical, resultsLowClinicalHighPRS = calculate_auc(df_binary, prs_methods, binary_measures)
#   results['auc_low_clinical'] = resultsLowClinical
#   results['auc_low_clinical_high_prs'] = resultsLowClinicalHighPRS
#   
#   # Calculate calibration
#   results_calibration = calculate_calibration(df_binary, prs_methods, binary_measures,'PHENOTYPE')
#   results['calibration'] = results_calibration
#   #results['auc_plots'] = auc_plots
#   
#   
#   results['nri'] = {}
#   for i, prs1 in enumerate(prs_methods):
#       results['nri'][prs1] = {}
#       for measure in binary_measures+['None']:
#           try:
#               prs_nri, prs_nri_cases, prs_nri_controls = calculate_single_prs_nri(df_binary,prs1,measure)
#               results['nri'][prs1][f"{measure}_low"] = {
#                   'nri': prs_nri,
#                   'nri_events': prs_nri_cases,
#                   'nri_non_events': prs_nri_controls
#               }
#           except Exception as e:
#               results['nri'][prs1][f'{measure}_low'] = {'error': str(e)}
#               
##           try:
##               prs_nri, prs_nri_cases, prs_nri_controls = calculate_single_prs_nri(df_binary,prs1,measure)
##               results['nri'][prs1][measure] = {
##                   'nri': prs_nri,
##                   'nri_events': prs_nri_cases,
##                   'nri_non_events': prs_nri_controls
##               }
##           except Exception as e:
##               results['nri'][prs1][measure] = {'error': str(e)}
#               
#           #calculate nri for prs compared to clinical measure
#           if measure != 'None':
#               try:
#                   results['nri'][f"{prs1}_vs_{measure}"] = {}
#                   nri, nri_events, nri_non_events = calculate_nri(
#                       df_binary, measure, prs1, "None", True
#                   )
#                   
#                   results['nri'][f"{prs1}_vs_{measure}"][measure] = {
#                       'nri': nri,
#                       'nri_events': nri_events,
#                       'nri_non_events': nri_non_events
#                   }
#                   
#               except Exception as e:
#                   results['nri'][f"{prs1}_vs_{measure}"][measure] = {'error': str(e)}
#                   
#               
#           for j, prs2 in enumerate(prs_methods):
#               if i >= j:  # Skip self-comparisons and repeated comparisons
#                   continue
#           
#               key = f"{prs2}_vs_{prs1}"
#               results['nri'][key] = {}
#   
#           
#               #get the calibration key in for comparing prs calculations
#               results['calibration'][key] = {}
#               
#   
#               try:
#                   nri, nri_events, nri_non_events = calculate_nri(
#                       df_binary, prs1, prs2, measure, False
#                   )
#                   results['nri'][key][f'{measure}_low'] = {
#                       'nri': nri,
#                       'nri_events': nri_events,
#                       'nri_non_events': nri_non_events
#                   }
#               except Exception as e:
#                   results['nri'][key][f'{measure}_low'] = {'error': str(e)}
#                   
#               try:    
#                   results['calibration'][key][f'{measure}_low'] = calculate_calibration_comparing_prs(df_binary[df_binary[measure] == 1],prs1,prs2)
#               except Exception as e:
#                   results['calibration'][key][f'{measure}_low'] = {'error': str(e)}
#
#   return results, df_binary_test, df_binary_validation

def compare_prs_performance(df, clinical_measures, figPath, file_ext, risk_thresholds=None,
                            outcome_column=None):                  
    
    """
    Comprehensive comparison of multiple PRS methods against multiple clinical measures
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing PRS scores and clinical outcome measures 

    clinical_measures : list
        List of column names for the different clinical outcome measures
    risk_thresholds : list or None
        Risk thresholds for NRI calculation. If None, quartiles of the first PRS method will be used.
    file_ext = str to concatenate to image file
        String of holdout or empty string if using validation set
    figPath = str
        file path for storage of heatmaps
    outcome_column : str or None
        Binary outcome column to use for determining risk direction
        
    Returns:
    --------
    dict
        Dictionary containing performance metrics
    pandas.DataFrame
        DataFrame with added binary columns
    """
    results = {}
    
    # First convert clinical measures to binary if needed
#   df_binary_test, df_binary_validation, used_thresholds, determined_directions = convert_to_binary(
#       df, df_validation, clinical_measures, thresholds=clinical_thresholds, 
#       high_risk_quintile=high_risk_quintile, risk_directions=risk_directions,
#       outcome_column=outcome_column
#   )
    
    df_copy = df.copy()
    # Create list of binary clinical measures
    binary_measures = [f"{measure}_binary" for measure in clinical_measures 
                       if f"{measure}_binary" in df_copy.columns]
    
    # If no binary measures were created, use original measures (assuming they're already binary)
    if not binary_measures:
        binary_measures = clinical_measures
        
        
    prs_methods = [col for col in df_copy.columns if 'prs' in col]
    
    
    # Calculate AUC for all methods and binary measures
    resultsLowClinical, resultsLowClinicalHighPRS = calculate_auc(df_copy, prs_methods, binary_measures,figPath,file_ext)
    results['auc_low_clinical'] = resultsLowClinical
    results['auc_low_clinical_high_prs'] = resultsLowClinicalHighPRS
    
    
    # Calculate calibration
#   results_calibration = calculate_calibration(df_binary, prs_methods, binary_measures,'PHENOTYPE')
#   results['calibration'] = results_calibration
    #results['auc_plots'] = auc_plots
    
    
    results['nri'] = {}
    for i, prs1 in enumerate(prs_methods):
        results['nri'][prs1] = {}
        for measure in binary_measures+['None']:
            
            #calculate nri for prs compared to clinical measure
            if measure != 'None':
                try:
                    results['nri'][f"{prs1}_vs_{measure}"] = {}
                    nri, nri_events, nri_non_events = calculate_nri(
                        df_copy, measure, prs1, "None", True
                    )
                    
                    results['nri'][f"{prs1}_vs_{measure}"][measure] = {
                        'nri': nri,
                        'nri_events': nri_events,
                        'nri_non_events': nri_non_events
                    }
                    df_copy[f'{prs1}_binary'] = (df[prs1] >= df[prs1].quantile(0.8)).astype(int)
                    fig = plot_reclassification_table(df_copy, measure, f'{prs1}_binary')
                    fig.savefig(f'{figPath}/reclassificationHeatMap.{measure}v{prs1}.{file_ext}.png',dpi=150, bbox_inches='tight')
                    plt.close(fig)
                    
                    fig = plot_low_clinical_reclassification_table(df_copy, measure, f'{prs1}_binary')
                    fig.savefig(f'{figPath}/reclassificationHeatMap.LowOnly_{measure}v{prs1}.{file_ext}.png',dpi=150, bbox_inches='tight')
                    plt.close(fig)
                    
                    
                except Exception as e:
                    results['nri'][f"{prs1}_vs_{measure}"][measure] = {'error': str(e)}
                    
            else: #if clinical measure == None, compare prs calculations to eachother    
                for j, prs2 in enumerate(prs_methods):
                    if i >= j:  # Skip self-comparisons and repeated comparisons
                        continue
            
                    key = f"{prs2}_vs_{prs1}"
                    results['nri'][key] = {}
            
            
                    #get the calibration key in for comparing prs calculations
#                   results['calibration'][key] = {}
            
            
                    try:
                        nri, nri_events, nri_non_events = calculate_nri(
                            df_copy, prs1, prs2, measure, False
                        )
                        results['nri'][key][f'{measure}_low'] = {
                            'nri': nri,
                            'nri_events': nri_events,
                            'nri_non_events': nri_non_events
                        }
                        df_copy[f'{prs1}_binary'] = (df[prs1] >= df[prs2].quantile(0.8)).astype(int)
                        df_copy[f'{prs2}_binary'] = (df[prs2] >= df[prs2].quantile(0.8)).astype(int)
                        fig = plot_reclassification_table(df_copy, f'{prs1}_binary', f'{prs2}_binary')
                        fig.savefig(f'{figPath}/reclassificationHeatMap.{prs1}v{prs2}.{file_ext}.png',dpi=150, bbox_inches='tight')
                        plt.close(fig)
                        
                    except Exception as e:
                        results['nri'][key][f'{measure}_low'] = {'error': str(e)}

                        
    return results



def impute_clinical_data(train_df, test_df, clinical_columns, 
                        prs_columns, outcome_column,
                        method='mean', visualize=True):
    """
    Impute missing values in clinical measures for PRS analysis
    
    Parameters:
    -----------
    train_df : pandas.DataFrame
        DataFrame containing clinical measures, PRS, and outcome variables for training data
    test_df : pandas.DataFrame
        DataFrame containing clinical measures, PRS, and outcome variables for training data
    clinical_columns : list
        List of clinical measure column names that might contain missing values
    prs_column : str
        Name of the PRS column
    outcome_column : str
        Name of the outcome column
    method : str
        Imputation method: 'iterative' (MICE), 'knn', 'mean', 'median', or 'regression'
    visualize : bool
        Whether to visualize the distribution before and after imputation
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with imputed values for clinical measures
    """

    # Make a copy of the original dataframe to avoid modifying it
    train_imputed = train_df.copy()
    test_imputed = test_df.copy()
    
    
    # Check for missing values
    missing_counts = train_df[clinical_columns].isnull().sum()
    missing_percent = 100 * missing_counts / len(train_df)
    
    print("Missing values in clinical measures:")
    for col, count, pct in zip(missing_counts.index, missing_counts, missing_percent):
        print(f"{col}: {count} values missing ({pct:.2f}%)")
        
    # Columns to use for imputation (can include PRS as predictor)
    predictor_cols = clinical_columns + prs_columns
    if outcome_column in train_df.columns:
        predictor_cols += [outcome_column]  # Including outcome can improve imputation
        
    # Extract only numeric columns for imputation
    numeric_cols = train_df[predictor_cols].select_dtypes(include=np.number).columns
        
    # Select imputation method
    if method == 'iterative':
        # MICE - Multivariate Imputation by Chained Equations
        imputer = IterativeImputer(
            estimator=RandomForestRegressor(n_estimators=100, random_state=42),
            random_state=42,
            max_iter=10
        )

            
        # Fit on training data
        imputer.fit(train_df[numeric_cols])
        
        for col in numeric_cols:
            if train_df[col].isnull().any():
                    train_mask = train_df[col].isnull()
                    train_imputed.loc[train_mask, col] = imputer.transform(train_df[numeric_cols])[train_mask, list(numeric_cols).index(col)]
                
            if test_df[col].isnull().any():
                    test_mask = test_df[col].isnull()
                    test_imputed.loc[test_mask, col] = imputer.transform(test_df[numeric_cols])[test_mask, list(numeric_cols).index(col)]
            
    elif method == 'knn':
        # K-Nearest Neighbors imputation
        imputer = KNNImputer(n_neighbors=5)
        
        # Fit on training data
        imputer.fit(train_df[numeric_cols])
        
        for col in numeric_cols:
            if train_df[col].isnull().any():
                    train_mask = train_df[col].isnull()
                    train_imputed.loc[train_mask, col] = imputer.transform(train_df[numeric_cols])[train_mask, list(numeric_cols).index(col)]
                
            if test_df[col].isnull().any():
                    test_mask = test_df[col].isnull()
                    test_imputed.loc[test_mask, col] = imputer.transform(test_df[numeric_cols])[test_mask, list(numeric_cols).index(col)]


            
    elif method in ['mean', 'median']:
        # Simple imputation with mean or median
        imputer = SimpleImputer(strategy=method)
        for col in clinical_columns:
            if train_df[col].isnull().any():
                imputer.fit(train_df[col].values.reshape(-1, 1) )
                
                train_imputed[col] = imputer.transform(train_df[col].values.reshape(-1, 1) )
                test_imputed[col] = imputer.transform(test_df[col].values.reshape(-1, 1) )
                
                
    elif method == 'regression':
        # Impute each variable using regression on other variables
        for target_col in clinical_columns:
            if train_df[target_col].isnull().any():
                # Identify rows with missing values in this column
                missing_mask = train_df[target_col].isnull()
                
                # Skip if all values are missing
                if missing_mask.all():
                    print(f"All values missing in {target_col}, using mean imputation")
                    train_imputed[target_col] = train_imputed[target_col].fillna(train_imputed[target_col].mean())

                    continue
                
                # Use other variables as predictors
                pred_cols = [c for c in predictor_cols if c != target_col]
                
                # Identify complete cases for training
                complete_cases = train_df[pred_cols + [target_col]].dropna().index
                
                # Skip if no complete cases
                if len(complete_cases) == 0:
                    print(f"No complete cases to train model for {target_col}, using mean imputation")
                    train_imputed[target_col] = train_df[target_col].fillna(train_df[target_col].mean())
                    
                    test_imputed[target_col] = test_df[target_col].fillna(train_df[target_col].mean())
                    continue
                
                # Train regression model
                model = RandomForestRegressor(n_estimators=100, random_state=42)
                model.fit(train_df.loc[complete_cases, pred_cols], train_df.loc[complete_cases, target_col])
                
                # Predict missing values
                missing_mask = test_df[target_col].isnull()
                missing_indices = test_df[missing_mask]
                
                predicted_values = model.predict(test_df.loc[missing_indices, pred_cols])
                test_imputed.loc[missing_indices, target_col] = predicted_values
                
    # Visualize the results if requested
    if visualize and any(missing_counts > 0):
        fig, axes = plt.subplots(len(clinical_columns), 2, figsize=(12, 4*len(clinical_columns)))
        if len(clinical_columns) == 1:
            axes = axes.reshape(1, 2)
            
        for i, col in enumerate(clinical_columns):
            if missing_counts[col] > 0:
                # Plot original distribution (excluding missing values)
                sns.histplot(test_df[col].dropna(), ax=axes[i, 0], kde=True)
                axes[i, 0].set_title(f"Original {col} Distribution (non-missing)")
                
                # Plot imputed distribution
                sns.histplot(test_imputed[col], ax=axes[i, 1], kde=True)
                axes[i, 1].set_title(f"Imputed {col} Distribution")
                
                # Compare with missing values highlighted
                missing_mask = test_df[col].isnull()
                if sum(missing_mask) > 0:
                    sns.scatterplot(x=test_imputed.index[missing_mask], 
                                     y=test_imputed.loc[missing_mask, col],
                                     color='red', ax=axes[i, 1], label='Imputed values')
                    axes[i, 1].legend()
            else:
                axes[i, 0].set_title(f"{col} - No missing values")
                axes[i, 1].set_title(f"{col} - No imputation needed")
                axes[i, 0].axis('off')
                axes[i, 1].axis('off')
                
        plt.tight_layout()
#       plt.show()
        
    # Calculate and print imputation quality metrics
    if any(missing_counts > 0):
        print("\nImputation summary:")
        for col in clinical_columns:
            if missing_counts[col] > 0:
                # Calculate basic statistics for the imputed values
                original_mean = test_df[col].mean()
                original_std = test_df[col].std()
                imputed_mean = test_imputed.loc[test_df[col].isnull(), col].mean()
                imputed_std = test_imputed.loc[test_df[col].isnull(), col].std()
                
                print(f"{col}:")
                print(f"  Original: mean={original_mean:.2f}, std={original_std:.2f}")
                print(f"  Imputed:  mean={imputed_mean:.2f}, std={imputed_std:.2f}")
                
    return train_imputed, test_imputed

def scale_data(df):

    # Initialize the StandardScaler
    scaler = StandardScaler()
    
    # Fit the scaler to your data (compute the mean and standard deviation)
    scaler.fit(df)
    
    # Transform the data using the fitted scaler
    scaled_data = scaler.transform(df)
    
    # Create a new DataFrame with the scaled data
    scaled_df = pd.DataFrame(scaled_data, columns=df.columns,index=df.index)

    return(scaled_df,scaler)

def calculate_idi(data, data_validation, outcome_column, prs_column, clinical_column):
    """
    Calculate Integrated Discrimination Improvement (IDI) for PRS added to clinical measures.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame containing the outcome variable, PRS, and clinical measures for training data
    scaled_data_validation: pandas.DataFrame
        DataFrame containing the outcome variable, PRS, and clinical measures for validation data
    outcome_column : str
        Name of the binary outcome column (0/1 or False/True)
    prs_column : str
        Name of the PRS column
    clinical_columns : list
        List of clinical measure column names to include in the reference model
    
    Returns:
    --------
    dict
        Dictionary containing IDI, its 95% CI, p-value, and component values
    """
    
    idi_results = {}
    
    if type(prs_column) == str:
        # Ensure all required columns exist
        required_columns = [clinical_column] + [outcome_column] + [prs_column]
    else:
        required_columns = [clinical_column] + [outcome_column] + prs_column
        
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"Required column '{col}' not found in data")
    
    scaled_data = data[data[f'{clinical_column}_binary'] == 1].copy()
    scaled_data_validation = data_validation[data_validation[f'{clinical_column}_binary'] == 1].copy()
    
    # Extract y (outcome)
    y = scaled_data[outcome_column].values
    
    # Create X matrices for both models
    X_clinical = scaled_data[[clinical_column]].values
    if type(prs_column) == str:
        X_combined = scaled_data[[clinical_column,prs_column]].values
    else:
        X_combined = scaled_data[[clinical_column]+prs_column].values
        
    
    # Fit logistic regression models
    clinical_model = LogisticRegression(random_state=42)
    combined_model = LogisticRegression(random_state=42)
    
    clinical_model.fit(X_clinical, y)
    combined_model.fit(X_combined, y)
    
    
    # Create X matrices for validation set
    X_clinical_validation = scaled_data_validation[[clinical_column]].values
    if type(prs_column) == str:
        X_combined_validation = scaled_data_validation[[clinical_column,prs_column]].values
    else:
        X_combined_validation = scaled_data_validation[[clinical_column]+prs_column].values

    y_validation = scaled_data_validation[outcome_column].values
    
    
    # Get predictions
    p_clinical = clinical_model.predict_proba(X_clinical_validation)[:, 1]
    p_combined = combined_model.predict_proba(X_combined_validation)[:, 1]
    
    
    # Separate cases and controls
    cases = y_validation == 1
    controls = y_validation == 0
    
    # Calculate discrimination slopes
    clinical_discr_slope = np.mean(p_clinical[cases]) - np.mean(p_clinical[controls])
    combined_discr_slope = np.mean(p_combined[cases]) - np.mean(p_combined[controls])
    
    # Calculate IDI
    idi = combined_discr_slope - clinical_discr_slope
    
    # Calculate components for interpretation
    idi_results[clinical_column] = {
        'IDI': idi,
        'clinical_auc': roc_auc_score(y_validation, p_clinical),
        'combined_auc': roc_auc_score(y_validation, p_combined),
        'mean_prob_cases_clinical': np.mean(p_clinical[cases]),
        'mean_prob_controls_clinical': np.mean(p_clinical[controls]),
        'mean_prob_cases_combined': np.mean(p_combined[cases]),
        'mean_prob_controls_combined': np.mean(p_combined[controls]),
        'clinical_discrimination_slope': clinical_discr_slope,
        'combined_discrimination_slope': combined_discr_slope
    }
    
    # Calculate confidence interval for IDI using bootstrap
    # (This is a simplified version - real applications should use more iterations)
    n_bootstrap = 1000
    idi_bootstrap = []
    
    for _ in range(n_bootstrap):
        # Bootstrap sample indices
        indices = np.random.choice(len(y_validation), len(y_validation), replace=True)
        
        # Calculate IDI on bootstrap sample
        y_boot = y_validation[indices]
        p_clinical_boot = p_clinical[indices]
        p_combined_boot = p_combined[indices]
        
        boot_cases = y_boot == 1
        boot_controls = y_boot == 0
        
        clinical_slope_boot = np.mean(p_clinical_boot[boot_cases]) - np.mean(p_clinical_boot[boot_controls])
        combined_slope_boot = np.mean(p_combined_boot[boot_cases]) - np.mean(p_combined_boot[boot_controls])
        
        idi_bootstrap.append(combined_slope_boot - clinical_slope_boot)
        
    # Calculate 95% confidence interval
    idi_bootstrap = np.array(idi_bootstrap)
    idi_results[clinical_column]['IDI_CI_lower'] = np.percentile(idi_bootstrap, 2.5)
    idi_results[clinical_column]['IDI_CI_upper'] = np.percentile(idi_bootstrap, 97.5)
    
    # Calculate p-value (proportion of bootstrap samples with IDI <= 0)
    idi_results[clinical_column]['p_value'] = np.mean(idi_bootstrap <= 0)
    
    fig = visualize_idi(p_clinical, p_combined, y_validation, idi_results[clinical_column],clinical_column,prs_column)
    
    return idi_results, fig

def visualize_idi(p_clinical, p_combined, y_validation,idi_results,clinical_marker,prs_method):
    """
    Visualize the IDI by plotting prediction probability distributions
    
    Parameters:
    -----------
    (Same as calculate_idi function)
    
    Returns:
    --------
    matplotlib.figure.Figure
        Figure containing the visualization
    """

    
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        'Outcome': y_validation,
        'Clinical Model': p_clinical,
        'Clinical + PRS Model': p_combined
    })
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot for clinical model
    sns.kdeplot(
        data=plot_df, x='Clinical Model', hue='Outcome', 
        palette={0:'blue', 1:'red'}, ax=axes[0], fill=True, alpha=0.5
    )
    axes[0].set_title('Clinical Model Predictions')
    axes[0].set_xlabel('Predicted Risk')
    axes[0].legend(['Cases', 'Controls'])
    
    # Plot for combined model
    sns.kdeplot(
        data=plot_df, x='Clinical + PRS Model', hue='Outcome',
        palette={0:'blue', 1:'red'}, ax=axes[1], fill=True, alpha=0.5
    )
    axes[1].set_title('Clinical + PRS Model Predictions')
    axes[1].set_xlabel('Predicted Risk')
    axes[1].legend(['Cases', 'Controls'])
    
    plt.tight_layout()
    
    # Calculate IDI for title
    fig.suptitle(f"IDI for {clinical_marker} using {prs_method} = {idi_results['IDI']:.4f} (95% CI: {idi_results['IDI_CI_lower']:.4f}-{idi_results['IDI_CI_upper']:.4f}, p={idi_results['p_value']:.4f})")
    plt.subplots_adjust(top=0.85)
    
    return fig

# Example usage
#def main(pheno,binary_to_use):
#   """Run an example using simulated data to demonstrate the function"""
#   # Generate simulated data
#
#
#   if binary_to_use:
#       file_ext = 'holdout'
#   else:
#       file_ext = 'validation'
#   
#   df = pd.read_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/reducedSHAP/combinedPRSGroupsWithcardio_main.csv')
#   
#   #prs columns to use in analysis
#   prs_columns = [col for col in df.columns if 'prs_' in col]
#   prs_columns = [col for col in prs_columns if 'scaled' not in col]
#
#   if pheno == 'type2Diabetes':
#       clinical_columns = ['Glycated haemoglobin (HbA1c)','Body mass index (BMI)','Glucose']
#       # Specify custom thresholds for some measures (optional)
#       clinical_thresholds = {
#           'Glycated haemoglobin (HbA1c)': 41,
#           'Body mass index (BMI)': 25,
#           'Glucose' : 12
#       }
#
#   else:
#       clinical_columns = ['Basal metabolic rate','Urea','Haemoglobin concentration']
#       clinical_thresholds = None
#       
#
#   
#   clinical_data = pd.read_csv('/Users/kerimulterer/ukbiobank/tanigawaData/environmentalCleaned/combinedmetabolomicsBloodChemistryFeatures.csv',usecols=['Participant ID']+clinical_columns)
#   clinical_data.rename(columns={'Participant ID':'IID'},inplace=True)
#   
#   data = df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
#   
#   
#   #download validation data to calculate 
#   validation_df = pd.read_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/reducedSHAP/holdout/CohortAssignedMax.EnvQuintiles.combinedPRSGroupsWithcardio_main1Kbins.holdout.csv')
#
#   
#   try:
#       validation_data = validation_df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
#   except Exception as e:
#       print('exception on merge : ',e)
#       prs_columns = [col for col in prs_columns if col in validation_df.columns]
#       validation_data = validation_df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
#       #get the column that wasn't in validation data
#       diff = list(set(list(data.columns)) - set(list(validation_data.columns)))
#       
#       data.drop(columns=diff,inplace=True)
#   
#   train_imputed, test_imputed = impute_clinical_data(data, validation_data, clinical_columns, prs_columns, 'PHENOTYPE', method='mean', visualize=False)
#   
#   scaled_data,scaler = scale_data(train_imputed.set_index(['IID','PHENOTYPE'])[prs_columns])
#   scaled_data.columns = prs_columns
#   scaled_data.reset_index(inplace=True)
#   
#   #merge with train_imputed 
#   scaled_data = scaled_data.merge(train_imputed[clinical_columns+['IID']], on=['IID'],how='left')
#   # Update values
#   scaled_data.loc[scaled_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
#   scaled_data.loc[scaled_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
#
#
#   
#   #scale validation data with scaler from training data
#   scaled_validation_data = scaler.transform(test_imputed.set_index(['IID','PHENOTYPE'])[prs_columns])
#   # Create a new DataFrame with the scaled data
#   scaled_validation_data = pd.DataFrame(scaled_validation_data, columns=prs_columns,index=validation_data.set_index(['IID','PHENOTYPE']).index)
#   scaled_validation_data.reset_index(inplace=True)
#   
#   #merge with train_imputed 
#   scaled_validation_data = scaled_validation_data.merge(test_imputed[clinical_columns + ['IID']], on=['IID'],how='left')
#   # Update values
#   scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
#   scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
#   
#   #merge the prscr_mix data not present in the validation data
#   if binary_to_use:
#       prs_validation = [col for col in validation_df.columns if 'scaled_prs' in col]
#       prs_test = [col for col in df.columns if 'scaled_prs' in col]
#       prs_diff = list(set(prs_validation)-set(prs_test))
#       scaled_validation_data = scaled_validation_data.merge(validation_df[['IID']+prs_diff],on=['IID'],how='left')
#       prs_columns = prs_columns + prs_diff
#   
#   
#   # Run comprehensive comparison with direction-aware threshold conversion
#   results, df_binary, df_validation_binary = compare_prs_performance(
#       scaled_data, scaled_validation_data, prs_columns, clinical_columns, 
#       binary_to_use=binary_to_use,
#       clinical_thresholds=clinical_thresholds,
#       outcome_column='PHENOTYPE'  # Use this to help determine risk direction
#   )
#   
#   if binary_to_use:
#       plot_df = df_validation_binary.copy()
#   else:
#       plot_df = df_binary.copy()
#
#   # Print thresholds and risk directions used
#   print("Thresholds and risk directions used for binarization:")
#   for measure in clinical_columns:
#       direction = "Higher values → higher risk" if results['risk_directions'][measure] else "Lower values → higher risk"
#       print(f"  {measure}: threshold = {results['thresholds_used'][measure]}, {direction}")
#       
#       # Plot risk distribution for clinical measures
#       fig = plot_risk_distribution(
#           plot_df, 
#           measure, 
#           f'{measure}_binary',
#           results['thresholds_used'][measure],
#           higher_is_riskier=results['risk_directions'][measure]
#       )
#       plt.title(f"Distribution of {measure}")
#       fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/riskDistribution.{measure}.{file_ext}.png')
#       
#       fig = plot_calibration_curves(results['calibration'], prs_columns, f'{measure}_binary', figsize=(10, 8))
#       fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/calibrationPlots.{measure}.lowRisk.{file_ext}.png')
#   
##       for prs in prs_columns:
##           if measure != 'None':
##               fig = results['auc_plots'][prs][f'{measure}_binary']
##           else:
##               fig = results['auc_plots'][prs][measure]
##           fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/AUC.{measure}.lowRisk.{prs}.png')
#           
#
#   
#   
#   # Print AUC results
#   print("\nAUC Results low clinical all prs:")
#   print(results['auc_low_clinical'])
##   results['auc_low_clinical'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossPRSLowClinicalMeasure.{file_ext}.csv')
#   
#   # Print AUC results
#   print("\nAUC Results low clinical high prs")
#   print(results['auc_low_clinical_high_prs'])
##   results['auc_low_clinical_high_prs'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossHighPRSLowClinicalMeasure.{file_ext}.csv')
#   
#   
#   # Print example NRI result
##   bin_measure = f"{measure}_binary"
##   print(f"\nNRI Example (Method1 vs Method2 for {bin_measure}):")
##   print(results['nri'][f'{prs_columns[0]}_vs_{prs_columns[1]}'][bin_measure])
#       
#
#   #calculate idi with both validation and holdout
#   if not binary_to_use: 
#   #   # Calculate IDI
#       results['idi'] = {}
#       
#       scaled_data,scaler = scale_data(train_imputed.set_index(['IID','PHENOTYPE']))
#       scaled_data.columns = prs_columns+clinical_columns
#       scaled_data.reset_index(inplace=True)
#       scaled_data.loc[scaled_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
#       scaled_data.loc[scaled_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
#   #   
#       #scale validation data with scaler from training data
#       scaled_validation_data = scaler.transform(test_imputed.set_index(['IID','PHENOTYPE']))
#       
#       # Create a new DataFrame with the scaled data
#       scaled_validation_data = pd.DataFrame(scaled_validation_data, columns=prs_columns+clinical_columns,index=validation_data.set_index(['IID','PHENOTYPE']).index)
#       scaled_validation_data.reset_index(inplace=True)
#       scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
#       scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
#       
#       for clinical_column in clinical_columns:
#           #attach the binary outcome for high, low risk
#           scaled_data = scaled_data.merge(df_binary[['IID',f'{clinical_column}_binary']],how='left')
#           scaled_validation_data = scaled_validation_data.merge(df_validation_binary[['IID',f'{clinical_column}_binary']],how='left')
#           
#           #calculate IDI for all prs calculations
#           idi_results,fig = calculate_idi(scaled_data, scaled_validation_data, 'PHENOTYPE', prs_columns, clinical_column)
#           results['idi']['all_prs_calculations'] = idi_results
#           
#           print("IDI Analysis Results: for combined PRS calculations and {clinical_column}")
#           print(f"IDI: {idi_results[clinical_column]['IDI']:.4f}")
#           print(f"95% CI: ({idi_results[clinical_column]['IDI_CI_lower']:.4f}, {idi_results[clinical_column]['IDI_CI_upper']:.4f})")
#           print(f"P-value: {idi_results[clinical_column]['p_value']:.4f}")
#           print(f"Clinical AUC: {idi_results[clinical_column]['clinical_auc']:.4f}")
#           print(f"Combined AUC: {idi_results[clinical_column]['combined_auc']:.4f}")
#           
#           fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/idi.{clinical_column}.performance.combinedPRS.{file_ext}.png')
#           
#           for prs_column in prs_columns:
#               idi_results,fig = calculate_idi(scaled_data, scaled_validation_data, 'PHENOTYPE', prs_column, clinical_column)
#               results['idi'][prs_column] = idi_results
#               
#               print(f"IDI Analysis Results: for prs {prs_column} and {clinical_column}")
#               print(f"IDI: {idi_results[clinical_column]['IDI']:.4f}")
#               print(f"95% CI: ({idi_results[clinical_column]['IDI_CI_lower']:.4f}, {idi_results[clinical_column]['IDI_CI_upper']:.4f})")
#               print(f"P-value: {idi_results[clinical_column]['p_value']:.4f}")
#               print(f"Clinical AUC: {idi_results[clinical_column]['clinical_auc']:.4f}")
#               print(f"Combined AUC: {idi_results[clinical_column]['combined_auc']:.4f}")
#               
#           
#               fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/idi.{clinical_column}.performance.{prs_column}.{file_ext}.png')
#
#
#   #put the auc results into a combined dataframe to save
#   clinical_auc = pd.DataFrame()
#   auc_list = [x for x in results.keys() if 'auc' in x]
#   for k in auc_list:
#       auc_df = results[k]
#       auc_df['threshold_used'] = k
#       clinical_auc = pd.concat([clinical_auc,auc_df],ignore_index=False)
#       del results[k]
#       
#   clinical_auc.reset_index(inplace=True)
#   clinical_auc.columns = ['prs_type'] + list(clinical_auc.columns[1:])
#   
#   clinical_auc.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/prsAUC.wrt.clinicalMeasures.{file_ext}.csv",index=False)
#   
#   # Save to a JSON file
##   with open(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/ClincalMeasuresValidationSet.json", "w") as f:
##       json.dump(results, f, indent=4)
#   
#   #put nri data into a csv
#   df_nri = pd.DataFrame()
#   for prs in results['nri'].keys():
#       prs_nri = pd.DataFrame(results['nri'][prs]).T
#       prs_nri['prs_calc'] = prs
#       prs_nri.reset_index(inplace=True)
#       prs_nri.columns = ['clinical_measure'] + list(prs_nri.columns[1:])
#       df_nri = pd.concat([df_nri,prs_nri],ignore_index=True)
#   df_nri.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/nri.wrt.clinicalMeasures.{file_ext}.csv",index=False)
#   
#
#   try:
#       #put idi values across clinical markers into a dataframe
#       df_idi = pd.DataFrame()
#       for prs in results['idi'].keys():
#           prs_idi = pd.DataFrame(results['idi'][prs]).T
#           prs_idi['prs_calc'] = prs
#           df_idi = pd.concat([df_idi,prs_idi],ignore_index=False)
#       df_idi.reset_index(inplace=True)
#       df_idi.columns = ['clinical_measure'] + list(df_idi.columns[1:])
#       df_idi.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/idi.wrt.clinicalMeasures.{file_ext}.csv",index=False)
#       
#   except KeyError:
#       print(f'no idi calculated for {binary_to_use} .... ')
#   
#   #get calibration data
#   df_calibration = pd.DataFrame()
#   for prs in results['calibration'].keys():
#       prs_cal = pd.DataFrame(results['calibration'][prs]).T
#       prs_cal['prs_calc'] = prs
#       df_calibration = pd.concat([df_calibration,prs_cal],ignore_index=False) 
#       
#   df_calibration.reset_index(inplace=True)
#   df_calibration.columns = ['clinical_measure'] + list(df_calibration.columns[1:])
#   df_calibration.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/calibration.wrt.clinicalMeasures.{file_ext}.csv",index=False)
#
#
#   for k in ['risk_directions','thresholds_used']:
#       try:
#           thresholds_risk[k] = pd.DataFrame(results[k],index=[k]).T
#       except UnboundLocalError:
#           thresholds_risk = pd.DataFrame(results[k],index=[k]).T
#           
#   thresholds_risk.reset_index(inplace=True)
#   thresholds_risk.columns = ['clinical_measure'] + list(thresholds_risk.columns[1:])
#           
#   thresholds_risk.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/thresholdsRiskDirection.wrt.clinicalMeasures.{file_ext}.csv",index=False)
    
    
#   return results, df_binary, df_validation_binary
    
def main(pheno,dataPath):
    """Run an example using simulated data to demonstrate the function"""
    
    
    
    # Generate simulated data
    figPath = f'{dataPath}/results/{pheno}/figures/clinicalFigures'
    scoresPath = f'{dataPath}/results/{pheno}/scores'
    
    df = pd.read_csv(f'{scoresPath}/combinedPRSGroups.csv')
    validation_df = pd.read_csv(f'{scoresPath}/combinedPRSGroups.holdout.csv')
    
    #prs columns to use in analysis
    prs_columns = [col for col in df.columns if 'prs_' in col]
    prs_columns = [col for col in prs_columns if 'scaled' not in col]
    
    if pheno == 'type2Diabetes':
        clinical_columns = ['Glycated haemoglobin (HbA1c)','Body mass index (BMI)','Glucose']
        # Specify custom thresholds for some measures (optional)
        clinical_thresholds = {
            'Glycated haemoglobin (HbA1c)': 41,
            'Body mass index (BMI)': 25,
            'Glucose' : 12
        }
        
    else:
        clinical_columns = ['Basal metabolic rate','Urea','Haemoglobin concentration']
        clinical_thresholds = None
        
        
    clinical_data = pd.read_csv(f'{dataPath}/participant_environment.csv',usecols=['Participant ID']+clinical_columns)
    clinical_data.rename(columns={'Participant ID':'IID'},inplace=True)
    
    data = df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
    
    
    try:
        validation_data = validation_df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
    except Exception as e:
        print('exception on merge : ',e)
        prs_columns = [col for col in prs_columns if col in validation_df.columns]
        validation_data = validation_df[['IID','PHENOTYPE']+prs_columns].merge(clinical_data, on=['IID'],how='left')
        #get the column that wasn't in validation data
        diff = list(set(list(data.columns)) - set(list(validation_data.columns)))
        
        data.drop(columns=diff,inplace=True)
        
    train_imputed, test_imputed = impute_clinical_data(data, validation_data, clinical_columns, prs_columns, 'PHENOTYPE', method='mean', visualize=False)
    
    scaled_data,scaler = scale_data(train_imputed.set_index(['IID','PHENOTYPE'])[prs_columns])
    scaled_data.columns = prs_columns
    scaled_data.reset_index(inplace=True)
    
    #merge with train_imputed 
    scaled_data = scaled_data.merge(train_imputed[clinical_columns+['IID']], on=['IID'],how='left')
    # Update values
    scaled_data.loc[scaled_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
    scaled_data.loc[scaled_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
    
    
    
    #scale validation data with scaler from training data
    scaled_validation_data = scaler.transform(test_imputed.set_index(['IID','PHENOTYPE'])[prs_columns])
    # Create a new DataFrame with the scaled data
    scaled_validation_data = pd.DataFrame(scaled_validation_data, columns=prs_columns,index=validation_data.set_index(['IID','PHENOTYPE']).index)
    scaled_validation_data.reset_index(inplace=True)
    
    #merge with train_imputed 
    scaled_validation_data = scaled_validation_data.merge(test_imputed[clinical_columns + ['IID']], on=['IID'],how='left')
    # Update values
    scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
    scaled_validation_data.loc[scaled_validation_data["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
    
    #merge the prscr_mix data not present in the validation data
    
    prs_validation = [col for col in validation_df.columns if 'scaled_prs' in col]
    prs_test = [col for col in df.columns if 'scaled_prs' in col]
    prs_diff = list(set(prs_validation)-set(prs_test))
    prs_diff = [col.replace('scaled_','') for col in prs_diff]
    scaled_validation_data = scaled_validation_data.merge(validation_df[['IID']+prs_diff],on=['IID'],how='left')
    
    # First convert clinical measures to binary if needed
    df_binary_test, df_binary_validation, used_thresholds, determined_directions = convert_to_binary(
        scaled_data, scaled_validation_data, clinical_columns, thresholds=clinical_thresholds, 
        high_risk_quintile=True, risk_directions=None,outcome_column='PHENOTYPE'
    )
    
    
    for item_tuple in [(df_binary_validation,'holdout'),(df_binary_test,'')]:
        
        
    # Print thresholds and risk directions used
        df_binary = item_tuple[0].copy()
        file_ext = item_tuple[1]
        
        # Run comprehensive comparison with direction-aware threshold conversion
        results = compare_prs_performance(
            df_binary, clinical_columns, figPath, file_ext,
            risk_thresholds=used_thresholds,
            outcome_column='PHENOTYPE'  # Use this to help determine risk direction
        )
        
        
        results['thresholds_used'] = used_thresholds
        results['risk_directions'] = determined_directions
        
        print("Thresholds used for binarization:")
        print(f'\nthresholds: {used_thresholds}')
        
        
        # Print AUC results
        print(f"\n Clinical measures for {file_ext} results")
        print("\nAUC Results low clinical all prs:")
        print(results['auc_low_clinical'])
#       results['auc_low_clinical'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossPRSLowClinicalMeasure.{file_ext}.csv')
        
        # Print AUC results
        print("\nAUC Results low clinical high prs")
        print(results['auc_low_clinical_high_prs'])
#       results['auc_low_clinical_high_prs'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossHighPRSLowClinicalMeasure.{file_ext}.csv')
        
        #put the auc results into a combined dataframe to save
        clinical_auc = pd.DataFrame()
        auc_list = [x for x in results.keys() if 'auc' in x]
        for k in auc_list:
            auc_df = results[k]
            auc_df['threshold_used'] = k
            clinical_auc = pd.concat([clinical_auc,auc_df],ignore_index=False)
            del results[k]
            
        clinical_auc.reset_index(inplace=True)
        clinical_auc.columns = ['prs_type'] + list(clinical_auc.columns[1:])
        
        clinical_auc.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/prsAUC.wrt.clinicalMeasures.{file_ext}.csv",index=False)
        
        for measure in clinical_columns:
            direction = "Higher values → higher risk" if results['risk_directions'][measure] else "Lower values → higher risk"
            print(f"  {measure}: threshold = {results['thresholds_used'][measure]}, {direction}")
            
            # Plot risk distribution for clinical measures
            fig = plot_risk_distribution(
                df_binary, 
                measure, 
                f'{measure}_binary',
                results['thresholds_used'][measure],
                higher_is_riskier=results['risk_directions'][measure]
            )
            plt.title(f"Distribution of {measure}")
            fig.savefig(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/clinicalFigures/riskDistribution.{measure}.{file_ext}.png')
            plt.close(fig)

            
        #put nri data into a csv
        df_nri = pd.DataFrame()
        for prs in results['nri'].keys():
            prs_nri = pd.DataFrame(results['nri'][prs]).T
            prs_nri['prs_calc'] = prs
            prs_nri.reset_index(inplace=True)
            prs_nri.columns = ['clinical_measure'] + list(prs_nri.columns[1:])
            df_nri = pd.concat([df_nri,prs_nri],ignore_index=True)
        df_nri.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/nri.wrt.clinicalMeasures.{file_ext}.csv",index=False)

        
        
        for k in ['risk_directions','thresholds_used']:
            try:
                thresholds_risk[k] = pd.DataFrame(results[k],index=[k]).T
            except UnboundLocalError:
                thresholds_risk = pd.DataFrame(results[k],index=[k]).T
                
        thresholds_risk.reset_index(inplace=True)
        thresholds_risk.columns = ['clinical_measure'] + list(thresholds_risk.columns[1:])
        
        thresholds_risk.to_csv(f"/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/data/clinicalMeasures/thresholdsRiskDirection.wrt.clinicalMeasures.{file_ext}.csv",index=False)
        
        
        
        
        
if __name__ == '__main__':
    
    pheno = 'type2Diabetes'
    main(pheno)
    
    

    
#if __name__ == '__main__':
#   
#   pheno = 'celiacDisease'
#   for t in [False,True]:
#       main(pheno,t)
    
    
    

    

    
    
    