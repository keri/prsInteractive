#!/usr/bin/env python3
"""
PRS Model Filtering Module

Identifies and filters redundant PRS models based on statistical metrics.
Uses multi-metric approach: Kappa, Jaccard, and DeLong test.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Set, Optional
import warnings
import re


def is_distinct_pair(kappa: float, jaccard: float, 
                     delong_pvalue: float, delta_auc: float,
                     kappa_threshold: float = 0.40,
                     jaccard_threshold: float = 0.50,
                     pvalue_threshold: float = 0.05,
                     auc_threshold: float = 0.05) -> Tuple[bool, str]:
    """
    Determine if two PRS models are statistically distinct using multi-metric criteria.
    
    Uses 2-out-of-3 rule:
    1. Low agreement (Kappa < threshold)
    2. Low overlap (Jaccard < threshold)  
    3. Different discrimination (DeLong p < threshold AND ΔAUC > threshold)
    
    Parameters:
    -----------
    kappa : float
        Cohen's Kappa coefficient
    jaccard : float
        Jaccard index (0-1, where 1 = complete overlap)
    delong_pvalue : float
        P-value from DeLong test
    delta_auc : float
        Absolute difference in AUC between models
    kappa_threshold : float, default=0.40
        Threshold for considering Kappa as indicating distinctness
    jaccard_threshold : float, default=0.50
        Threshold for considering Jaccard as indicating distinctness
    pvalue_threshold : float, default=0.05
        P-value threshold for DeLong test significance
    auc_threshold : float, default=0.05
        Minimum AUC difference to consider meaningful
        
    Returns:
    --------
    Tuple[bool, str]
        (is_distinct, reason_string)
        - is_distinct: True if models should both be kept, False if redundant
        - reason_string: Explanation of the decision
    """
    criteria_met = []
    criteria_failed = []
    
    # Criterion 1: Low agreement (Kappa)
    if pd.isna(kappa) or kappa == 'NA':
        kappa_val = None
    else:
        kappa_val = float(kappa)
        
    if kappa_val is not None and kappa_val < kappa_threshold:
        criteria_met.append(f"Kappa={kappa_val:.3f} < {kappa_threshold} (low agreement)")
    elif kappa_val is not None:
        criteria_failed.append(f"Kappa={kappa_val:.3f} ≥ {kappa_threshold} (moderate/high agreement)")
    
    # Criterion 2: Low overlap (Jaccard)
    if pd.isna(jaccard) or jaccard == 'NA':
        jaccard_val = None
    else:
        jaccard_val = float(jaccard)
        
    if jaccard_val is not None and jaccard_val < jaccard_threshold:
        criteria_met.append(f"Jaccard={jaccard_val:.3f} < {jaccard_threshold} ({100*(1-jaccard_val):.1f}% unique)")
    elif jaccard_val is not None:
        criteria_failed.append(f"Jaccard={jaccard_val:.3f} ≥ {jaccard_threshold} (high overlap)")
    
    # Criterion 3: Different discrimination (DeLong + ΔAUC)
    if pd.isna(delong_pvalue) or delong_pvalue == 'NA':
        delong_pval = None
    else:
        delong_pval = float(delong_pvalue)
        
    if pd.isna(delta_auc) or delta_auc == 'NA':
        delta_auc_val = None
    else:
        delta_auc_val = abs(float(delta_auc))
    
    if (delong_pval is not None and delta_auc_val is not None and 
        delong_pval < pvalue_threshold and delta_auc_val > auc_threshold):
        criteria_met.append(
            f"DeLong p={delong_pval:.4f}, delta AUC={delta_auc_val:.3f} (significantly different)"
        )
    elif delong_pval is not None and delta_auc_val is not None:
        criteria_failed.append(
            f"DeLong p={delong_pval:.4f}, delta AUC={delta_auc_val:.3f} (not significantly different)"
        )
    
    # Decision: Keep both if ≥2 criteria met
    n_criteria_met = len(criteria_met)
    is_distinct = n_criteria_met >= 2
    
    if is_distinct:
        reason = f"DISTINCT ({n_criteria_met}/3 criteria): " + "; ".join(criteria_met)
    else:
        reason = f"REDUNDANT ({n_criteria_met}/3 criteria): " + "; ".join(criteria_failed)
    
    return is_distinct, reason

def _resolve_redundant_pair(model: str, redundant_model: str, model_performance: dict) -> tuple[str, str, str]:
    """
    Decide which of two redundant models to filter out.

    Priority:
      1. Substring rule  – if one model's name is wholly contained within the
         other's, the longer (composite) name is filtered.  This removes
         'epi+main' when it is redundant with 'epi' or 'main', and generalises
         to any future composite naming scheme.
      2. Performance fallback – keep whichever has the higher AUC.

    Returns
    -------
    (to_filter, to_keep, reason_str)
    """
    # --- Rule 1: substring containment ---
    if redundant_model in model:          # e.g. 'epi' in 'epi+main'
        return model, redundant_model, f"substring rule ('{redundant_model}' is a component of '{model}')"
    if model in redundant_model:          # e.g. 'epi+main' contains 'epi'
        return redundant_model, model, f"substring rule ('{model}' is a component of '{redundant_model}')"
        
    # --- Rule 2: performance fallback ---
    if model_performance[model] >= model_performance[redundant_model]:
        return redundant_model, model, "performance fallback (higher AUC kept)"
    else:
        return model, redundant_model, "performance fallback (higher AUC kept)"


def filter_redundant_models(stats_results: pd.DataFrame,
                            model_column: str = 'comparison',
                            performance_metric: str = 'auc',
                            verbose: bool = True) -> Dict:
    """
    Identify redundant PRS models from statistical comparison results.

    Analyzes all pairwise comparisons and returns the minimal set of distinct
    models that should be kept, filtering out redundant ones.

    Parameters
    ----------
    stats_results : pd.DataFrame
        DataFrame with statistical test results. Must contain columns:
        - comparison : e.g. 'main v epi'
        - test       : name of test; rows where test is 'Cohen_Kappa',
                       'Jaccard_Overlap', or 'DeLong' are used
        - statistic  : test statistic value
        - pvalue     : p-value (for DeLong rows)
        - auc1, auc2 : AUC values (for DeLong rows)
    model_column : str, default 'comparison'
        Column containing model comparison strings (e.g. 'main v epi').
    performance_metric : str, default 'auc'
        Metric used for ranking when the performance fallback is triggered.
    verbose : bool, default True
        Print detailed filtering decisions.

    Returns
    -------
    Dict with keys:
        'models_to_keep'    : List of model names to retain
        'models_to_filter'  : List of model names to remove (redundant)
        'pairwise_decisions': DataFrame with all pairwise distinctness decisions
        'redundancy_groups' : Dict mapping filtered models to what they are
                              redundant with, plus the resolution reason
        'model_performance' : Dict of model -> max AUC seen across comparisons
    """
  
    # ------------------------------------------------------------------ #
    # 1. Extract all unique model names from comparison strings            #
    # ------------------------------------------------------------------ #
    all_models = set()
    for comp in stats_results[model_column].unique():
        if ' v ' in comp:
            all_models.update(m.strip() for m in comp.split(' v '))
          
    all_models = sorted(all_models)
  
    if verbose:
        print(f"\n{'='*70}")
        print("PRS Model Redundancy Analysis")
        print(f"{'='*70}")
        print(f"Found {len(all_models)} models: {', '.join(all_models)}")
        print(f"{'='*70}\n")
      
    # ------------------------------------------------------------------ #
    # 2. Build pairwise comparison table                                   #
    # ------------------------------------------------------------------ #
    pairwise_results = []
  
    for comp in stats_results[model_column].unique():
        if ' v ' not in comp:
            continue
      
        parts = comp.split(' v ')
        if len(parts) != 2:
            continue
      
        model1, model2 = [m.strip() for m in parts]
        comp_data = stats_results[stats_results[model_column] == comp]
      
        # Cohen's Kappa
        kappa_row = comp_data[comp_data['test'] == 'Cohen_Kappa']
        kappa = kappa_row['statistic'].values[0] if len(kappa_row) > 0 else np.nan
      
        # Jaccard Overlap
        jaccard_row = comp_data[comp_data['test'] == 'Jaccard_Overlap']
        jaccard = jaccard_row['statistic'].values[0] if len(jaccard_row) > 0 else np.nan
      
        # DeLong
        delong_row = comp_data[comp_data['test'] == 'DeLong']
        if len(delong_row) > 0:
            delong_p = delong_row['pvalue'].values[0]
            auc1 = delong_row['auc1'].values[0] if 'auc1' in delong_row.columns else np.nan
            auc2 = delong_row['auc2'].values[0] if 'auc2' in delong_row.columns else np.nan
            delta_auc = abs(auc1 - auc2) if not (pd.isna(auc1) or pd.isna(auc2)) else np.nan
        else:
            delong_p = auc1 = auc2 = delta_auc = np.nan
          
        is_distinct, reason = is_distinct_pair(kappa, jaccard, delong_p, delta_auc)
      
        pairwise_results.append({
            'model1': model1,
            'model2': model2,
            'comparison': comp,
            'kappa': kappa,
            'jaccard': jaccard,
            'delong_p': delong_p,
            'auc1': auc1,
            'auc2': auc2,
            'delta_auc': delta_auc,
            'is_distinct': is_distinct,
            'reason': reason,
        })
      
        if verbose:
            status = "✅ KEEP BOTH" if is_distinct else "❌ REDUNDANT"
            print(f"{model1} vs {model2}: {status}")
            print(f"  {reason}\n")
          
    pairwise_df = pd.DataFrame(pairwise_results)
  
    # ------------------------------------------------------------------ #
    # 3. Build redundancy graph                                            #
    # ------------------------------------------------------------------ #
    redundancy_graph: Dict[str, list] = {}
  
    for _, row in pairwise_df[~pairwise_df['is_distinct']].iterrows():
        m1, m2 = row['model1'], row['model2']
        redundancy_graph.setdefault(m1, []).append(m2)
        redundancy_graph.setdefault(m2, []).append(m1)
      
    # ------------------------------------------------------------------ #
    # 4. Collect per-model performance (max AUC seen across comparisons)  #
    # ------------------------------------------------------------------ #
    model_performance: Dict[str, float] = {}
    for model in all_models:
        aucs = []
        for _, row in pairwise_df.iterrows():
            if row['model1'] == model and not pd.isna(row['auc1']):
                aucs.append(row['auc1'])
            if row['model2'] == model and not pd.isna(row['auc2']):
                aucs.append(row['auc2'])
        model_performance[model] = max(aucs) if aucs else 0.0
      
    # ------------------------------------------------------------------ #
    # 5. Resolve redundancies                                              #
    #    Priority:                                                         #
    #      (a) Substring rule  – composite/longer name is filtered first  #
    #      (b) Performance fallback – lower AUC is filtered               #
    # ------------------------------------------------------------------ #
    models_to_keep = set(all_models)
    models_to_filter: Dict[str, Dict] = {}   # filtered_model -> {kept, reason}
  
    for model, redundant_peers in redundancy_graph.items():
        if model not in models_to_keep:
            continue  # already resolved
      
        for redundant_model in redundant_peers:
            if redundant_model not in models_to_keep:
                continue  # already resolved
          
            to_filter, to_keep, res_reason = _resolve_redundant_pair(
                model, redundant_model, model_performance
            )
          
            models_to_keep.discard(to_filter)
            models_to_filter[to_filter] = {'kept': to_keep, 'reason': res_reason}
          
            if to_filter == model:
                break  # this model is now filtered; skip remaining peers
          
    models_to_keep_sorted = sorted(models_to_keep)
  
    # ------------------------------------------------------------------ #
    # 6. Verbose summary                                                   #
    # ------------------------------------------------------------------ #
    if verbose:
        print(f"\n{'='*70}")
        print("FILTERING SUMMARY")
        print(f"{'='*70}")
        print(f"\n✅ MODELS TO KEEP ({len(models_to_keep_sorted)}):")
        for model in models_to_keep_sorted:
            auc = model_performance.get(model, 'N/A')
            auc_str = f"{auc:.4f}" if isinstance(auc, float) else auc
            print(f"  • {model} (AUC: {auc_str})")
          
        if models_to_filter:
            print(f"\n❌ MODELS TO FILTER ({len(models_to_filter)}):")
            for filtered, info in models_to_filter.items():
                kept = info['kept']
                res_reason = info['reason']
                auc_f = model_performance.get(filtered, 'N/A')
                auc_k = model_performance.get(kept, 'N/A')
                auc_f_str = f"{auc_f:.4f}" if isinstance(auc_f, float) else auc_f
                auc_k_str = f"{auc_k:.4f}" if isinstance(auc_k, float) else auc_k
                print(f"  • {filtered} (AUC: {auc_f_str}) → redundant with {kept} "
                      f"(AUC: {auc_k_str}) [{res_reason}]")
        else:
            print("\n✅ All models are distinct – no filtering needed!")
          
        print(f"{'='*70}\n")
      
    return {
        'models_to_keep': models_to_keep_sorted,
        'models_to_filter': list(models_to_filter.keys()),
        'pairwise_decisions': pairwise_df,
        'redundancy_groups': models_to_filter,   # now includes 'reason' per entry
        'model_performance': model_performance,
    }


def apply_model_filter(
  data: pd.DataFrame, 
  models_to_keep: List[str],
  model_columns: Optional[List[str]] = None,
  base_columns: Optional[List[str]] = None,
  verbose: bool = True
) -> pd.DataFrame:
  """
  Filter a dataset to only include specified PRS models.
  
  Correctly handles compound model names (e.g., 'epi+main') to avoid
  false matches and properly excludes columns from models not in models_to_keep.
  
  Parameters:
  -----------
  data : pd.DataFrame
    Original dataset with PRS scores
  models_to_keep : List[str]
    List of model names to retain (e.g., ['cardio', 'epi', 'main'])
  model_columns : List[str], optional
    Explicit list of model column names to keep.
    If None, auto-detects using pattern matching.
  base_columns : List[str], optional
    List of base column names to always keep (e.g., ['IID', 'PHENOTYPE']).
    If None, auto-detects common base columns.
  verbose : bool, default=True
    Print filtering information
    
  Returns:
  --------
  pd.DataFrame
    Filtered dataset with only kept models
    
  Examples:
  ---------
  >>> models = ['cardio', 'epi', 'main']  # Without 'epi+main'
  >>> filtered_df = apply_model_filter(df, models)
  # Keeps: scaled_prs_cardio, bin_epi, etc.
  # Excludes: scaled_prs_epi+main, bin_epi+main, etc.
  """
  # Define base columns to always keep
  if base_columns is None:
    base_columns = ['IID', 'PHENOTYPE']
    
  # Keep base columns that exist
  base_cols_present = [col for col in base_columns if col in data.columns]
  
  # Sort models by length (longest first) to avoid partial matches
  models_sorted = sorted(models_to_keep, key=len, reverse=True)
  
  if model_columns is None:
    # Common PRS column patterns
    # These identify ANY model column, not just ones we want to keep
    model_patterns = [
      r'scaled_prs_',
      r'prs_',
      r'bin_',
      r'percentile_',
      r'score_',
      r'risk_',
    ]
    
    # Auto-detect model columns
    model_columns_to_keep = []
    
    for col in data.columns:
      # Skip base columns
      if col in base_cols_present:
        continue
      
      # Check if this is a model column (matches any pattern)
      is_model_column = any(re.search(pattern, col) for pattern in model_patterns)
      
      if not is_model_column:
        continue  # Not a model column, skip it
      
      # This IS a model column - check if it matches a model we want to keep
      matched = False
      for model in models_sorted:
        # Pattern: model name with delimiters (_, -, or start/end)
        # Matches: 'scaled_prs_epi', 'bin_main'
        # Does NOT match: 'epi' in 'scaled_prs_epi+main'
        pattern = rf'(^|_|-){re.escape(model)}(_|-|$)'
        
        if re.search(pattern, col):
          model_columns_to_keep.append(col)
          matched = True
          break  # Found match, stop checking
        
      # If it's a model column but doesn't match any model to keep,
      # it will be excluded (not added to model_columns_to_keep)
        
    model_columns = model_columns_to_keep
    
  # Identify true non-model columns (columns that aren't base columns and aren't model columns)
  model_patterns = [
    r'scaled_prs_',
    r'prs_',
    r'bin_',
    r'percentile_',
    r'score_',
    r'risk_',
  ]
  
  other_cols = []
  for col in data.columns:
    # Skip if it's a base column
    if col in base_cols_present:
      continue
    
    # Skip if already identified as a model column to keep
    if col in model_columns:
      continue
    
    # Check if it's ANY model column (to be excluded)
    is_any_model_column = any(re.search(pattern, col) for pattern in model_patterns)
    
    if not is_any_model_column:
      # This is a non-model column (metadata, etc.)
      other_cols.append(col)
      
  # Combine: base columns + model columns we want + other non-model columns
  cols_to_keep = base_cols_present + model_columns + other_cols
  
  # Remove duplicates while preserving order
  cols_to_keep = list(dict.fromkeys(cols_to_keep))
  
  # Filter dataframe
  filtered_data = data[cols_to_keep].copy()
  
  if verbose:
    print(f"\n{'='*70}")
    print(f"Model Column Filtering")
    print(f"{'='*70}")
    print(f"Models to keep: {', '.join(models_to_keep)}")
    print(f"Original columns: {len(data.columns)}")
    print(f"Filtered columns: {len(filtered_data.columns)}")
    print(f"Removed columns: {len(data.columns) - len(filtered_data.columns)}")
    
    # Show columns by model
    print(f"\nModel columns kept:")
    for model in models_to_keep:
      pattern = rf'(^|_|-){re.escape(model)}(_|-|$)'
      model_cols = [col for col in model_columns if re.search(pattern, col)]
      if model_cols:
        print(f"  {model}: {len(model_cols)} columns")
        for col in sorted(model_cols)[:5]:
          print(f"    - {col}")
        if len(model_cols) > 5:
          print(f"    ... and {len(model_cols) - 5} more")
          
    # Show base columns
    if base_cols_present:
      print(f"\nBase columns: {', '.join(base_cols_present)}")
      
    # Show other columns
    if other_cols:
      print(f"\nOther columns: {len(other_cols)}")
      if len(other_cols) <= 5:
        print(f"  {', '.join(other_cols)}")
      else:
        print(f"  {', '.join(other_cols[:5])}, ... ({len(other_cols) - 5} more)")
        
    # Show removed columns (model columns not in models_to_keep)
    removed_cols = [col for col in data.columns if col not in cols_to_keep]
    if removed_cols:
      print(f"\nRemoved model columns: {len(removed_cols)}")
      if len(removed_cols) <= 10:
        print(f"  {', '.join(removed_cols)}")
      else:
        print(f"  {', '.join(removed_cols[:5])}, ... ({len(removed_cols) - 5} more)")
        
    print(f"{'='*70}\n")
    
  return filtered_data

# Convenience function for common use case
def quick_filter_models(stats_csv_path: str, 
                       exclude_comparisons: List[str] = None,
                       output_path: str = None) -> Dict:
    """
    Quick filtering from CSV file path.
    
    Parameters:
    -----------
    stats_csv_path : str
        Path to CSV with statistical test results
    exclude_comparisons : List[str], optional
        List of comparison patterns to exclude (e.g., ['all v combined'])
    output_path : str, optional
        Path to save filtering results CSV
        
    Returns:
    --------
    Dict with filtering results
    """
    # Load data
    stats = pd.read_csv(stats_csv_path)
    
    # Filter out excluded comparisons
    if exclude_comparisons:
        for pattern in exclude_comparisons:
            stats = stats[~stats['comparison'].str.contains(pattern, case=False, na=False)]
    
    # Run filtering
    results = filter_redundant_models(stats, verbose=True)
    
    # Save results if requested
    if output_path:
        results['pairwise_decisions'].to_csv(output_path, index=False)
        print(f"Saved pairwise decisions to: {output_path}")
    
    return results


if __name__ == "__main__":
    # Example usage
    print("PRS Model Filtering Module")
    print("="*70)
    print("\nExample usage:")
    print("""
    from prs_model_filter import filter_redundant_models, quick_filter_models
    
    # Quick filtering from CSV
    results = quick_filter_models(
        'McNemarStatsTestsAcrossPRSCalculations_refactored.csv',
        exclude_comparisons=['all v combined'],
        output_path='pairwise_filtering_decisions.csv'
    )
    
    # Get models to keep
    models_to_keep = results['models_to_keep']
    print(f"Keep these models: {models_to_keep}")
    
    # Apply filter to your PRS data
    from prs_model_filter import apply_model_filter
    
    prs_data = pd.read_csv('your_prs_scores.csv')
    filtered_data = apply_model_filter(prs_data, models_to_keep)
    filtered_data.to_csv('filtered_prs_scores.csv', index=False)
    """)