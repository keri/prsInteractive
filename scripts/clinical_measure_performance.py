#!/usr/bin/env python3

#Keri Multerer April 7, 2025
#methods to calculate PRS performance against clinical measures

import numpy as np
import json
import os
import argparse
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
from helper.draw_plots import *


COLOR_SCHEME = {
    'case-control': {
        'cases': '#dfa6a5',      # Purple/Pink
        'controls': '#b7e1f3',      # Blue
        'total': '#916F6F',   # Teal/Green
        'name': 'case-control'
    }}

# Extended color scheme for product/summed variants
COHORT_COLORS_EXTENDED = {
    'main': '#E69F00',      # Orange
    'epi': '#56B4E9',       # Sky blue
    'epi+main': '#CC79A7',  # Pinkish purple
    'cardio': '#009E73',    # Bluish green
    'all': '#F0E442',       # Yellow
    'combined': '#D55E00',  # Vermillion
    # Product variants (same as base)
    'epi_product': '#56B4E9',
    'cardio_product': '#009E73',
    'epi+main_product': '#CC79A7',
    'all_product': '#F0E442',
    'main_product': '#E69F00',
    # Summed variants (darker shades)
    'epi_summed': '#0073A8',      # Darker blue
    'cardio_summed': '#006B52',   # Darker green
    'epi+main_summed': '#9F5580', # Darker purple
    'all_summed': '#C7B800',      # Darker yellow
    'main_summed': '#C77D00',     # Darker orange
}



def extract_prs_models_from_data(df, exclude_combined=True):
    """
    Dynamically extract PRS model names from dataframe columns.
    
    Looks for *_high_risk columns to identify available PRS models.
    Returns base model names without suffixes (e.g., 'main', 'epi_product', 'cardio_summed').
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing PRS columns
    exclude_combined : bool, default=True
        Whether to exclude 'combined' from the list
        
    Returns
    -------
    list
        List of PRS model names found in the data
    """
    prs_models = []
    
    # Look for _high_risk columns
    high_risk_cols = [col for col in df.columns if col.endswith('_high_risk')]
    
    for col in high_risk_cols:
        model_name = col.replace('_high_risk', '')
        # Exclude combined, any, and all models
        if exclude_combined and 'combined' in model_name.lower():
            continue
        if model_name.lower() in ['any', 'all']:
            continue
        prs_models.append(model_name)
    
    # If no _high_risk columns found, try extracting from scaled_prs_ or bin_ columns
    if not prs_models:
        for col in df.columns:
            if col.startswith('scaled_prs_'):
                model_name = col.replace('scaled_prs_', '')
                if 'threshold' not in model_name:
                    if exclude_combined and 'combined' in model_name.lower():
                        continue
                    prs_models.append(model_name)
            elif col.startswith('bin_'):
                model_name = col.replace('bin_', '')
                if exclude_combined and 'combined' in model_name.lower():
                    continue
                prs_models.append(model_name)
    
    # Remove duplicates and sort
    prs_models = sorted(set(prs_models))
    
    return prs_models


def extract_base_cohorts_from_models(prs_models):
    """
    Extract cohort names from PRS model names, keeping product/summed as separate cohorts.
    
    For stacked bars, we want each variant as its own cohort (e.g., 'epi_product' and 
    'epi_summed' stay separate, not grouped under 'epi').
    
    Parameters
    ----------
    prs_models : list
        List of PRS model names (e.g., ['main', 'epi_product', 'cardio_summed'])
        
    Returns
    -------
    dict
        Dictionary with full model names as keys, empty lists as values
        e.g., {'main': [], 'epi_product': [], 'epi_summed': [], 'cardio_product': [], 'cardio_summed': []}
    """
    cohorts = []
    
    for model in prs_models:
        # Skip compound models with '+' 
        if '+' in model:
            continue
        
        # Keep the full model name (including _product/_summed suffixes)
        cohorts.append(model)
    
    # Create dictionary with empty lists
    return {cohort: [] for cohort in sorted(cohorts)}


def create_sankey_plot_clinical_data(df, figPath, use_epi_main=False):
    
    # Dynamically extract PRS models from the data
    prs_methods = extract_prs_models_from_data(df, exclude_combined=True)
    
    # Legacy parameter for backward compatibility
    if not use_epi_main:
        # Filter out models containing 'epi+main' or '+' compounds if use_epi_main is False
        prs_methods = [m for m in prs_methods if '+' not in m]
    
    if not prs_methods:
        print("WARNING: No PRS models found in data")
        return
    
    print(f"Using PRS methods: {prs_methods}")
    print(f"Total holdout samples: {len(df)}")
    print(f"Cases in holdout: {df['PHENOTYPE'].sum()}")
    
    # check which binary columns exist:
    binary_cols = [col for col in df.columns if col.endswith('_binary')]
    print(f"\nAvailable binary columns: {binary_cols}")
    
    
    for binar_col in binary_cols:
        
        # Filter for cases with low clinical risk (binary = 0)
        df_filtered = df[(df['PHENOTYPE'] == 1) & (df[binar_col] == 0)].copy()
        
        # Track which individuals have already been counted
        counted_individuals = set()
        
        # Count flows from each PRS method to combined
        flows = []
        for method in prs_methods:
            high_risk_col = f'{method}_high_risk'
            
            if high_risk_col not in df_filtered.columns:
                print(f"Warning: {high_risk_col} not found in data")
                continue
            
            # Get individuals high-risk in this method who haven't been counted yet
            high_risk_mask = (df_filtered[high_risk_col] == 1)
            uncounted_mask = ~df_filtered.index.isin(counted_individuals)
            current_method_mask = high_risk_mask & uncounted_mask
            
            # Count individuals high-risk in this method (not yet counted) -> combined high-risk
            high_to_high = len(df_filtered[current_method_mask & 
                                            (df_filtered['combined_high_risk'] == 1)])
            
#           # Count individuals high-risk in this method -> combined low-risk
#           high_to_low = len(df_filtered[current_method_mask & 
#                                           (df_filtered['combined_high_risk'] == 0)])
                                        
            # Add these individuals to the counted set
            counted_individuals.update(df_filtered[current_method_mask].index)
            
            # Add flows only if there are individuals
            if high_to_high > 0:
                flows.append({
                    'source': f'{method.upper()} High Risk',
                    'target': 'Combined High Risk',
                    'value': high_to_high,
                    'color': COHORT_COLORS[method]
                })
            
#           if high_to_low > 0:
#               flows.append({
#                   'source': f'{method.upper()} High Risk',
#                   'target': 'Combined Low Risk',
#                   'value': high_to_low,
#                   'color': COHORT_COLORS[method]
#               })
                
            print(f"{method.upper()}: {df_filtered[high_risk_col].sum()} high-risk individuals")
            print(f"  -> Combined High: {high_to_high}")
#           print(f"  -> Combined Low: {high_to_low}")
            
            print(f"\nCombined high-risk: {df_filtered['combined_high_risk'].sum()}")
        
        print(f"\nTotal individuals counted: {len(counted_individuals)}")
        print(f"Combined high-risk total: {df_filtered['combined_high_risk'].sum()}")
            
        # Create Sankey diagram
        # Build unique labels
        all_labels = []
        label_dict = {}
        
        for flow in flows:
            if flow['source'] not in label_dict:
                label_dict[flow['source']] = len(all_labels)
                all_labels.append(flow['source'])
            if flow['target'] not in label_dict:
                label_dict[flow['target']] = len(all_labels)
                all_labels.append(flow['target'])
                
        # Build source, target, value, and color lists
        sources = [label_dict[flow['source']] for flow in flows]
        targets = [label_dict[flow['target']] for flow in flows]
        values = [flow['value'] for flow in flows]
        link_colors = [flow['color'] for flow in flows]
        
        # Create node colors
        node_colors = []
        node_x = []  # x position (0 = left, 1 = right)
        node_y = []  # y position (0 = top, 1 = bottom)
        
        for i, label in enumerate(all_labels):
            # Set x position
            if 'Combined' in label:
                node_x.append(1.0)  # Right side
            else:
                node_x.append(0.0)  # Left side
                
            # Set y position and color
            if 'CARDIO' in label:
                node_colors.append(COHORT_COLORS_EXTENDED['cardio'])
                node_y.append(0.1)  # Top
            elif 'EPI' in label and 'main' not in label.lower():
                node_colors.append(COHORT_COLORS_EXTENDED['epi'])
                node_y.append(0.4)  # Middle
            elif 'MAIN' in label:
                node_colors.append(COHORT_COLORS_EXTENDED['main'])
                node_y.append(0.7)  # Bottom
            elif 'Combined Low' in label:
                node_colors.append('#999999')  # Gray for combined low risk
                node_y.append(0.2)  # Top on right side
            elif 'Combined High' in label:
                node_colors.append('#D55E00')  # Red-orange for combined high risk
                node_y.append(0.6)  # Bottom on right side
            else:
                node_colors.append('#CCCCCC')
                node_y.append(0.5)
                    
                    
        # Create Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color='black', width=0.5),
#               label=all_labels,
                color=node_colors,
                x=node_x,
                y=node_y
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
            )
        )])
        
        fig.update_layout(
            title=dict(
                text=f"PRS High-Risk Classification Flow to Combined PRS<br>" + 
                    f"<sub>Cases with Low Clinical Risk (n={len(df_filtered)})</sub>",
                x=0.5,
                xanchor='center'
            ),
            font=dict(size=12, family='Arial'),
            width=1000,
            height=600,
            margin=dict(l=20, r=20, t=80, b=20)
        )
        
        # Save figure
        fig.write_html(f'{figPath}/Sankey{binar_col}.PRStoCombinedHighRisk.html')
        print(f"\nSankey plot saved as 'Sankey{binar_col}.PRStoCombinedHighRisk.html'")
        
        pio.write_image(fig,f'{figPath}/Sankey{binar_col}.PRStoCombinedHighRisk.png')
        
        # Also display if in interactive environment
        fig.show()
        
        # Print summary statistics
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        
        for method in prs_methods:
            high_risk_col = f'{method}_high_risk'
            if high_risk_col in df_filtered.columns:
                n_high = df_filtered[high_risk_col].sum()
                overlap_combined = len(df_filtered[(df_filtered[high_risk_col] == 1) & 
                                                    (df_filtered['combined_high_risk'] == 1)])
                if n_high > 0:
                    overlap_pct = (overlap_combined / n_high) * 100
                    print(f"{method.upper()}: {overlap_pct:.1f}% of high-risk -> Combined high-risk")
                    
        combined_high = df_filtered['combined_high_risk'].sum()
        print(f"\nTotal combined high-risk: {combined_high} ({combined_high/len(df_filtered)*100:.1f}%)")
        
        
        
def plot_nri_from_reclassification(df, clinical_vars, prs_col, show_total_as_line=False):
    
    """
    Create NRI bar chart showing correct reclassification for cases and controls
    across multiple clinical variables.
    
    Clinical measures: 0=Low Risk, 1=High Risk
    PRS: 0=Low Risk, 1=High Risk
    
    When prs_col contains 'combined', creates stacked bars showing the contribution
    of each cohort (main, epi, cardio) to NRI for cases and controls separately.
    
    For each clinical variable:
    - Controls: Count those moved from clinical high risk (1) to PRS low risk (0)
        (i.e., clinical says high risk [1], PRS correctly says low risk [0])
    - Cases: Count those moved from clinical low risk (0) to PRS high risk (1)
        (i.e., clinical says low risk [0], PRS correctly says high risk [1])
    
    Parameters:
    -----------
    df : pandas.DataFrame
            DataFrame with 'PHENOTYPE' column (1=case, 2=control)
    clinical_vars : list of clinical var str
            list clinical risk category column names (0=low risk, 1=high risk)
    prs_col : str
            Column name for PRS risk categories (0=low risk, 1=high risk)
    show_total_as_line : bool, default=False
            Ignored when prs_col contains 'combined'
            
    Returns:
    --------
    fig : matplotlib.figure.Figure
            The generated figure
    nri_summary : pandas.DataFrame
            Summary table with NRI values and counts
    """
    # Check if this is a combined PRS analysis
    is_combined = 'combined' in prs_col.lower()
    
    # Store results for each clinical variable
    results = []
    
    # Split by phenotype
    cases = df['PHENOTYPE'] == 1
    controls = df['PHENOTYPE'] == 0
    
    for risk_cat1 in clinical_vars:
        
        # =====================================================================
        # CONTROLS: Clinical says HIGH risk (1), want PRS to correctly say LOW risk (0)
        # =====================================================================
        reclass_controls = pd.crosstab(
            df.loc[controls, risk_cat1],  # Clinical risk (rows): 0=low, 1=high
            df.loc[controls, prs_col],     # PRS risk (columns): 0=low, 1=high
            rownames=[risk_cat1],
            colnames=[prs_col]
        )
        
        # Count controls correctly reclassified: 
        # Clinical HIGH risk (1) → PRS LOW risk (0)
        if 1 in reclass_controls.index and 0 in reclass_controls.columns:
            controls_correct = reclass_controls.loc[1, 0]
        else:
            controls_correct = 0
            
        # Total controls at high risk in clinical measure (clinical = 1)
        total_controls_high = reclass_controls.loc[1].sum() if 1 in reclass_controls.index else 0
        
        # NRI- for controls (proportion correctly kept at/moved to low risk)
        nri_controls = controls_correct / len(df[controls]) if len(df[controls]) > 0 else 0
        
        # =====================================================================
        # CASES: Clinical says LOW risk (0), want PRS to correctly say HIGH risk (1)
        # =====================================================================
        reclass_cases = pd.crosstab(
            df.loc[cases, risk_cat1],      # Clinical risk (rows): 0=low, 1=high
            df.loc[cases, prs_col],         # PRS risk (columns): 0=low, 1=high
            rownames=[risk_cat1],
            colnames=[prs_col]
        )
        
        # Count cases correctly reclassified:
        # Clinical LOW risk (0) → PRS HIGH risk (1)
        if 0 in reclass_cases.index and 1 in reclass_cases.columns:
            cases_correct = reclass_cases.loc[0, 1]
        else:
            cases_correct = 0
            
        # Total cases at low risk in clinical measure (clinical = 0)
        total_cases_low = reclass_cases.loc[0].sum() if 0 in reclass_cases.index else 0
        
        # NRI+ for cases (proportion correctly moved to high risk)
        nri_cases = cases_correct / len(df[cases]) if len(df[cases]) > 0 else 0
        
        # =====================================================================
        # COHORT-SPECIFIC NRI (only for combined PRS)
        # =====================================================================
        cohort_nri = {}
        if is_combined:
            # ── Dynamically discover available PRS models - keep as separate cohorts ──
            # Extract all available PRS models from the data
            available_models = extract_prs_models_from_data(df, exclude_combined=True)
            
            # Build cohort groups - each model is its own cohort now
            base_groups = extract_base_cohorts_from_models(available_models)
            
            # Populate groups with their corresponding _high_risk columns using EXACT matching
            for col in df.columns:
                if not col.endswith('_high_risk') or 'combined' in col:
                    continue
                model_name = col.replace('_high_risk', '')   # e.g. 'epi_product'
                
                # Exact match only (no prefix matching)
                if model_name in base_groups:
                    base_groups[model_name].append(col)
            
            # Define priority order for stacking (highest to lowest priority)
            # Each case is counted only once, attributed to highest-priority cohort
            priority_order = ['main', 'epi_product', 'epi_summed', 'cardio_product', 'cardio_summed']
            
            # Sort cohorts by priority (models not in priority_order go at the end alphabetically)
            sorted_cohorts = []
            for p in priority_order:
                if p in base_groups:
                    sorted_cohorts.append(p)
            # Add any remaining cohorts not in priority list
            for cohort in sorted(base_groups.keys()):
                if cohort not in sorted_cohorts:
                    sorted_cohorts.append(cohort)

            def _any_high(group_cols):
                if not group_cols:
                    return pd.Series(False, index=df.index)
                return df[group_cols].max(axis=1).astype(bool)

            # Create cohort-specific high-risk masks using priority order
            cohort_highs = {}
            for cohort_name in sorted_cohorts:
                cohort_highs[cohort_name] = _any_high(base_groups[cohort_name])

            # For CASES: Clinical LOW (0) → PRS HIGH (1)
            cases_reclass_mask = (df[risk_cat1] == 0) & (df[prs_col] == 1) & cases

            # Apply priority: first cohort in priority order gets highest priority
            # Build exclusion masks for priority-based attribution (ensures each case counted ONCE)
            cases_by_cohort = {}
            exclusion_mask = pd.Series(False, index=df.index)
            
            for cohort_name in sorted_cohorts:
                cohort_high = cohort_highs[cohort_name]
                cases_by_cohort[cohort_name] = (cohort_high & cases_reclass_mask & ~exclusion_mask).sum()
                exclusion_mask = exclusion_mask | cohort_high

            # For CONTROLS: proportional attribution by cohort size
            controls_reclass_mask  = (df[risk_cat1] == 1) & (df[prs_col] == 0) & controls
            controls_reclass_total = controls_reclass_mask.sum()

            # Calculate total for each cohort (using priority order for consistency)
            totals_by_cohort = {cohort: cohort_highs[cohort].sum() 
                               for cohort in sorted_cohorts}
            total_all = sum(totals_by_cohort.values())

            # Proportional attribution for controls
            controls_by_cohort = {}
            if total_all > 0:
                for cohort_name in sorted_cohorts:
                    controls_by_cohort[cohort_name] = (
                        controls_reclass_total * (totals_by_cohort[cohort_name] / total_all)
                    )
            else:
                controls_by_cohort = {cohort: 0 for cohort in sorted_cohorts}

            n_cases    = len(df[cases])
            n_controls = len(df[controls])

            # Build cohort_nri dynamically using priority order
            cohort_nri = {}
            for cohort_name in sorted_cohorts:
                cohort_nri[f'cases_{cohort_name}'] = (
                    cases_by_cohort[cohort_name] / n_cases if n_cases > 0 else 0
                )
                cohort_nri[f'cases_{cohort_name}_n'] = cases_by_cohort[cohort_name]
                cohort_nri[f'controls_{cohort_name}'] = (
                    controls_by_cohort[cohort_name] / n_controls if n_controls > 0 else 0
                )
                cohort_nri[f'controls_{cohort_name}_n'] = int(controls_by_cohort[cohort_name])
            
        # =====================================================================
        # Store results
        # =====================================================================
        result_dict = {
            'variable': risk_cat1.replace('_risk_cat', '').replace('_', ' '),
            'nri_cases': nri_cases,
            'nri_controls': nri_controls,
            'total_nri': nri_cases + nri_controls,
            'cases_correct_n': cases_correct,
            'cases_total_low_n': total_cases_low,
            'controls_correct_n': controls_correct,
            'controls_total_high_n': total_controls_high,
            'total_cases': len(df[cases]),
            'total_controls': len(df[controls])
        }
        
        if is_combined:
            result_dict.update(cohort_nri)
            
        results.append(result_dict)
        
    # Convert to DataFrame
    nri_summary = pd.DataFrame(results)
    
    # =========================================================================
    # CREATE THE BAR CHART
    # =========================================================================
    
    # Get selected color scheme
    colors = COLOR_SCHEME['case-control']
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Set positions for bars
    x = np.arange(len(nri_summary))
    
    if is_combined:
        # =====================================================================
        # STACKED BARS: Cases and Controls NRI by cohort (DYNAMIC)
        # =====================================================================
        width = 0.35
        
        # Detect which cohorts are in nri_summary and preserve priority order
        cohort_cols = [col for col in nri_summary.columns if col.startswith('cases_') and not col.endswith('_n')]
        all_cohorts = [col.replace('cases_', '') for col in cohort_cols]
        
        # Apply priority ordering to cohorts
        priority_order = ['main', 'epi_product', 'epi_summed', 'cardio_product', 'cardio_summed']
        cohorts = []
        for p in priority_order:
            if p in all_cohorts:
                cohorts.append(p)
        # Add any remaining cohorts not in priority list (alphabetically)
        for c in sorted(all_cohorts):
            if c not in cohorts:
                cohorts.append(c)
        
        if not cohorts:
            print("WARNING: No cohort-specific data found for combined PRS")
            # Fall back to non-stacked layout
            is_combined = False
        else:
            print(f"Creating stacked bars for cohorts: {cohorts}")
            
            # CASES STACKED BAR (left) - build dynamically
            cases_bars = {}
            bottom_cases = pd.Series(0, index=nri_summary.index)
            
            for cohort in cohorts:
                cases_col = f'cases_{cohort}'
                if cases_col in nri_summary.columns:
                    # Get color for this cohort (fallback to gray if not found)
                    color = COHORT_COLORS_EXTENDED.get(cohort, '#CCCCCC')
                    
                    bars = ax.bar(
                        x - width/2,
                        nri_summary[cases_col],
                        width,
                        bottom=bottom_cases,
                        label=f'{cohort.capitalize()} PRS' if cohort == cohorts[0] else None,
                        color=color,
                        alpha=0.85,
                        edgecolor='black',
                        linewidth=0.5
                    )
                    cases_bars[cohort] = bars
                    bottom_cases = bottom_cases + nri_summary[cases_col]
            
            # CONTROLS STACKED BAR (right) - build dynamically
            controls_bars = {}
            bottom_controls = pd.Series(0, index=nri_summary.index)
            
            for cohort in cohorts:
                controls_col = f'controls_{cohort}'
                if controls_col in nri_summary.columns:
                    color = COHORT_COLORS_EXTENDED.get(cohort, '#CCCCCC')
                    
                    bars = ax.bar(
                        x + width/2,
                        nri_summary[controls_col],
                        width,
                        bottom=bottom_controls,
                        color=color,
                        alpha=0.85,
                        edgecolor='black',
                        linewidth=0.5
                    )
                    controls_bars[cohort] = bars
                    bottom_controls = bottom_controls + nri_summary[controls_col]
            
            # Add total NRI labels above each bar
            for i in range(len(x)):
                # Cases total
                cases_total = nri_summary.iloc[i]['nri_cases']
                ax.text(
                    x[i] - width/2, cases_total + 0.01,
                    f'NRI+: {cases_total:.3f}',
                    ha='center', va='bottom',
                    fontsize=9, fontweight='bold',
                    color=colors['cases']
                )
                
                # Controls total
                controls_total = nri_summary.iloc[i]['nri_controls']
                ax.text(
                    x[i] + width/2, controls_total + 0.01,
                    f'NRI-: {controls_total:.3f}',
                    ha='center', va='bottom',
                    fontsize=9, fontweight='bold',
                    color=colors['controls']
                )
    
    if not is_combined:
        # =====================================================================
        # ORIGINAL NON-STACKED LAYOUT
        # =====================================================================
        width = 0.25
        
        # Plot NRI+ for cases
        bars_cases = ax.bar(
            x - width,
            nri_summary['nri_cases'],
            width,
            label=f'NRI+ (Cases: Clinical Low→PRS High risk)\nn={nri_summary["total_cases"].iloc[0]:,}',
            color=colors['cases'],
            alpha=0.8,
            linewidth=1
        )
        
        # Plot NRI- for controls
        bars_controls = ax.bar(
            x,
            nri_summary['nri_controls'],
            width,
            label=f'NRI- (Controls: Clinical High→PRS Low risk)\nn={nri_summary["total_controls"].iloc[0]:,}',
            color=colors['controls'],
            alpha=0.8,
            linewidth=1
        )
        
        # Plot total NRI
        if show_total_as_line:
            ax.plot(
                x + width,
                nri_summary['total_nri'],
                marker='o',
                markersize=8,
                linewidth=2,
                color=colors['total'],
                label='Total NRI',
                zorder=5
            )
            for i, (xi, val) in enumerate(zip(x + width, nri_summary['total_nri'])):
                ax.text(
                    xi, val + 0.01,
                    f'{val:.3f}',
                    ha='center', va='bottom',
                    fontsize=10, fontweight='bold',
                    color='#1f77b4'
                )
        else:
            bars_total = ax.bar(
                x + width,
                nri_summary['total_nri'],
                width,
                label='Total NRI',
                color=colors['total'],
                alpha=0.8,
                linewidth=1
            )
            # Add value labels
            for bar in bars_total:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{height:.3f}',
                    ha='center', va='bottom',
                    fontsize=10, fontweight='bold'
                )
                
        # Add value labels on NRI bars
        for i, bar in enumerate(bars_cases):
            height = bar.get_height()
            n_correct = nri_summary.iloc[i]['cases_correct_n']
            n_total = nri_summary.iloc[i]['cases_total_low_n']
            ax.text(
                bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}\n({n_correct}/{n_total})',
                ha='center', va='bottom',
                fontsize=8
            )
            
        for i, bar in enumerate(bars_controls):
            height = bar.get_height()
            n_correct = nri_summary.iloc[i]['controls_correct_n']
            n_total = nri_summary.iloc[i]['controls_total_high_n']
            ax.text(
                bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}\n({n_correct}/{n_total})',
                ha='center', va='bottom',
                fontsize=8
            )
            
    # Customize the plot
    ax.set_xlabel('Clinical Risk Variables', fontsize=13, fontweight='bold')
    ax.set_ylabel('Net Reclassification Improvement (NRI)', fontsize=13, fontweight='bold')
    
    if is_combined:
        title_text = (f'Cohort Contribution to NRI: {prs_col.replace("_", " ")}\n'
                     'Stacked bars show Main→Epi→Cardio (priority) contribution to correct reclassification\n'
                     'Clinical: 1=Low Risk, 0=High Risk | PRS: 1=High Risk, 0=Low Risk')
    else:
        title_text = (f'Net Reclassification Improvement: Clinical Risk vs {prs_col.replace("_", " ")}\n'
                     'Clinical: 1=Low Risk, 0=High Risk | PRS: 1=High Risk, 0=Low Risk')
        
    ax.set_title(title_text, fontsize=13, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(nri_summary['variable'], fontsize=11)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Legend
    if is_combined:
        ax.legend(loc='upper left', fontsize=10, frameon=True, shadow=True)
    else:
        ax.legend(loc='upper left', fontsize=10, frameon=True, shadow=True)
        
    # Set y-axis limits with padding
    y_max = max(nri_summary['nri_cases'].max(), nri_summary['nri_controls'].max())
    ax.set_ylim(-0.02, y_max * 1.25)
    
    plt.tight_layout()
    
    return fig, nri_summary

    
def plot_reclassification_table(df, risk_cat1, risk_cat2,show_total_as_line=False):
    """
    Create side-by-side reclassification heatmaps for cases and controls.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with 'PHENOTYPE' column (1=case, 0=control)
    risk_cat1 : str (will be clinical measure)
        Column name for method 1 risk categories
    risk_cat2 : str  
        Column name for method 2 risk categories

show_total_as_line : bool, default=False
    If True, show total NRI as a line plot instead of bars
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # =========================================================================
    # CASES Panel
    # =========================================================================
    cases = df['PHENOTYPE'] == 1
    

    reclass_cases = pd.crosstab(
        df.loc[cases, risk_cat1],  
        df.loc[cases, risk_cat2],   
        rownames=[risk_cat1],
        colnames=[risk_cat2],
        margins=True
    )
    
    sns.heatmap(
        reclass_cases, 
        annot=True, fmt='d', cmap='Blues', 
        ax=axes[0], cbar=False,
        linewidths=0.5, linecolor='gray'
    )
    axes[0].set_title('Reclassification Table: CASES', fontsize=14, fontweight='bold')
    axes[0].set_xlabel(f'{risk_cat2}\n(0=Low Risk, 1=High Risk)', fontsize=11)  
#   if 'prs' in risk_cat1:
    axes[0].set_ylabel(f'{risk_cat1}\n(0=Low Risk, 1=High Risk)', fontsize=11)  
#   else:
#       axes[0].set_ylabel(f'{risk_cat1}\n(1=Low Risk, 0=High Risk)', fontsize=11) 
        
        
        
    # =========================================================================
    # CONTROLS Panel
    # =========================================================================
    controls = df['PHENOTYPE'] == 0
    
    reclass_controls = pd.crosstab(
        df.loc[controls, risk_cat1], 
        df.loc[controls, risk_cat2],  
        rownames=[risk_cat1],         
        colnames=[risk_cat2],         
        margins=True
    )
    
    
    sns.heatmap(
        reclass_controls, 
        annot=True, fmt='d', cmap='Greens',
        ax=axes[1], cbar=False,
        linewidths=0.5, linecolor='gray'
    )
    axes[1].set_title('Reclassification Table: CONTROLS', fontsize=14, fontweight='bold')
    axes[1].set_xlabel(f'{risk_cat2}\n(0=Low Risk, 1=High Risk)', fontsize=11)  
#   if 'prs' in risk_cat1:
    axes[0].set_ylabel(f'{risk_cat1}\n(0=Low Risk, 1=High Risk)', fontsize=11)  
#   else:
#       axes[0].set_ylabel(f'{risk_cat1}\n(1=Low Risk, 0=High Risk)', fontsize=11)  
        
    plt.tight_layout()
    
    return fig

def plot_low_clinical_reclassification_table(df, risk_cat1, risk_cat2):
    """
    Create reclassification heatmap for cases and controls within low risk clinical.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with 'PHENOTYPE' column (1=case, 0=control)
    risk_cat1 : str
        Column name for method 1 risk categories
    risk_cat2 : str  
        Column name for method 2 risk categories
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # =========================================================================
    # CASES Panel
    # =========================================================================
    
    df_low = df[df[risk_cat1] == 0]
    
    controls = df_low['PHENOTYPE'] == 0
    cases = df_low['PHENOTYPE'] == 1
    
    
    # ✅ Fast indexing - no chained indexing
    reclass_low_risk_clinical = pd.crosstab(
        df_low['PHENOTYPE'],   
        df_low[risk_cat2],   
        rownames=['PHENOTYPE'],
        colnames=[risk_cat2],
        margins=True
    )
    
    sns.heatmap(
        reclass_low_risk_clinical, 
        annot=True, fmt='d', cmap='Greys', 
        cbar=False,
        linewidths=0.5, linecolor='gray'
    )
    ax.set_title('Reclassification Table: LOW CLINICAL RISK', fontsize=14, fontweight='bold')
    ax.set_xlabel(f'{risk_cat2}\n(0=Low Risk, 1=High Risk)', fontsize=11)  # ✅ cat2 on x-axis
    ax.set_ylabel(f'T2D in low {risk_cat1}\n(0=Controls, 1=Cases)', fontsize=11)  # ✅ cat2 on x-axis
    
    
    
    plt.tight_layout()
    return fig



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
    
    # Set quintile threshold (default to 0.8 or 0.2 if True), else use thresholds set with known clinical measures used in risk assessment
    if high_risk_quintile is True:
        quintile_threshold = 0.2
    elif isinstance(high_risk_quintile, (int, float)) and 1 <= high_risk_quintile <= 0:
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
                df_result[f"{measure}_binary"] = (df_result[measure] >= threshold).astype(int)
                df_result_validation[f"{measure}_binary"] = (df_result_validation[measure] >= threshold).astype(int)
            else:
                df_result[f"{measure}_binary"] = (df_result[measure] <= threshold).astype(int)
                df_result_validation[f"{measure}_binary"] = (df_result_validation[measure] <= threshold).astype(int)
                
            used_thresholds[measure] = threshold
        else:
            print(f"Warning: No threshold specified for '{measure}' and quintile option not used")
            
    return df_result, df_result_validation, used_thresholds, determined_directions


def calculate_auc(dfFull, prs_method, prs_binary, clinical_measures):
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
    resultsAllPRS = pd.DataFrame(index=[prs_method], columns=clinical_measures)
    resultsTop20PRS = pd.DataFrame(index=[prs_method], columns=clinical_measures)
    

#   for prs in prs_methods:
    for measure in clinical_measures:
        #get the low risk clinical measure
        df = dfFull[dfFull[measure] == 0]
        
        #get the high prs 
#       threshold = df[prs].quantile(.80)
#       df['high_prs'] = (df[prs_method] <= threshold).astype(int)
        
        # Ensure the clinical measure is binary (0/1)
        if not set(df['PHENOTYPE'].unique()).issubset({0, 1}):
            raise ValueError(f"PHENOTYPE must be binary (0/1)")
            
        # Calculate AUC
        try:
            auc = roc_auc_score(df['PHENOTYPE'], df[prs_method])
            resultsAllPRS.loc[prs_method, measure] = auc
        except Exception as e:
            resultsAllPRS.loc[prs_method, measure] = np.nan
            print(f"Error calculating AUC for {prs_method} vs low risk {measure}: {e}")
            
#       try:
#           if 'scaled_prs' in prs_method:
#               threshold = df[prs_method].quantile(.80)
#               df['high_prs'] = (df[prs_method] >= threshold).astype(int)
#           else:
#               df['high_prs'] = (df[prs_method] >= 800).astype(int)
        dfHigh = df[df[prs_binary] == 1]
        auc = roc_auc_score(dfHigh['PHENOTYPE'], dfHigh[prs_method])
        try:
            resultsTop20PRS.loc[prs_method, measure] = auc
            
        except Exception as e:
            resultsTop20PRS.loc[prs_method, measure] = np.nan
            print(f"Error calculating AUC for high risk {prs_method} vs low risk {measure}: {e}")
                
    return resultsAllPRS, resultsTop20PRS


def _count_reclassifications(old_risk: pd.Series, new_risk: pd.Series,
                             is_case: pd.Series) -> dict:
    """
    Core reclassification counter shared by calculate_nri and
    net_reclassification_index.

    Both caller functions resolve their inputs to binary (0/1 or bool) series
    before calling this helper, keeping counting logic in one place.

    Parameters
    ----------
    old_risk : pd.Series
        Binary (0/1 or bool) high-risk classification under the reference model.
    new_risk : pd.Series
        Binary (0/1 or bool) high-risk classification under the new model.
        Must share the same index as old_risk.
    is_case : pd.Series
        Boolean mask; True for cases (events), False for controls (non-events).
        Must share the same index as old_risk.

    Returns
    -------
    dict with keys:
        up_cases, down_cases       - raw counts for events
        up_controls, down_controls - raw counts for non-events
        n_cases, n_controls        - group sizes
    """
    is_control = ~is_case

    # Events: moving to high risk is good (cases should be flagged)
    up_cases   = int(np.sum(~old_risk[is_case]  &  new_risk[is_case]))   # Low->High (good)
    down_cases = int(np.sum( old_risk[is_case]  & ~new_risk[is_case]))   # High->Low (bad)

    # Non-events: moving to low risk is good (controls should not be flagged)
    up_controls   = int(np.sum(~old_risk[is_control] &  new_risk[is_control]))  # Low->High (bad)
    down_controls = int(np.sum( old_risk[is_control] & ~new_risk[is_control]))  # High->Low (good)

    return {
        'up_cases':      up_cases,
        'down_cases':    down_cases,
        'up_controls':   up_controls,
        'down_controls': down_controls,
        'n_cases':       int(is_case.sum()),
        'n_controls':    int(is_control.sum()),
    }


def calculate_nri(df, prs_method1, prs_method2, clinical=True):
    """
    Calculate Net Reclassification Improvement (NRI) comparing two binary risk
    classifications.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'PHENOTYPE' (1=case, 0=control) and both method columns.
    prs_method1 : str
        Column name for the reference classification (binary 0/1).
    prs_method2 : str
        Column name for the new classification (binary 0/1).
    clinical : bool
        When True, validates that prs_method1 is binary before proceeding.

    Returns
    -------
    tuple : (nri, nri_events, nri_non_events)
    """
    if clinical:
        if not set(df[prs_method1].dropna().unique()).issubset({0, 1}):
            raise ValueError(
                f"Column '{prs_method1}' must be binary (0/1) when clinical=True. "
                f"Found values: {sorted(df[prs_method1].dropna().unique())}"
            )

    counts = _count_reclassifications(
        old_risk=df[prs_method1].astype(bool),
        new_risk=df[prs_method2].astype(bool),
        is_case=df['PHENOTYPE'] == 1,
    )

    n_cases    = counts['n_cases']
    n_controls = counts['n_controls']

    # NRI components expressed as proportions (Pencina 2008)
    nri_events     = (counts['up_cases']      - counts['down_cases'])    / n_cases
    nri_non_events = (counts['down_controls'] - counts['up_controls'])   / n_controls
    nri            = nri_events + nri_non_events

    return nri, nri_events, nri_non_events

def net_reclassification_index(data: pd.DataFrame, prs1_name: str,
                               prs2_name: str, threshold: int,
                               phenotype_col: str = 'PHENOTYPE') -> dict:
    """
    Calculate Net Reclassification Index (NRI) and Integrated Discrimination
    Improvement (IDI) with full statistics (SE, CI, p-value).

    Delegates reclassification counting to _count_reclassifications(), the same
    helper used by calculate_nri(), so both functions are guaranteed to apply
    identical counting logic.

    Parameters
    ----------
    data : pd.DataFrame
        Must contain '{phenotype_col}' (1=case, 0=control), 'bin_{prs1_name}',
        'bin_{prs2_name}', 'scaled_prs_{prs1_name}', 'scaled_prs_{prs2_name}'.
    prs1_name : str
        Base name of the reference PRS model (e.g. 'main').
    prs2_name : str
        Base name of the new PRS model to evaluate.
    threshold : int
        Inclusive lower bound for high-risk classification on the bin_ columns.
        Individuals with bin >= threshold are classed as high risk.
    phenotype_col : str
        Column name for the outcome; expected coding is 1=case, 0=control.

    Returns
    -------
    dict with keys: test, statistic, pvalue, conf_int, nri_events,
                    nri_non_events, idi
    """
    print(f"\n{'='*70}")
    print("Net Reclassification Index (NRI)")
    print(f"{'='*70}")

    is_case = data[phenotype_col] == 1

    # Resolve bin_ columns to boolean high-risk flags using inclusive >= threshold.
    # Using >= ensures the boundary bin is not silently dropped into the low-risk
    # group (old code used strict >, which excluded the boundary value).
    old_risk = (data[f'bin_{prs1_name}'] >= threshold)
    new_risk = (data[f'bin_{prs2_name}'] >= threshold)

    counts = _count_reclassifications(old_risk, new_risk, is_case)

    up_cases      = counts['up_cases']
    down_cases    = counts['down_cases']
    up_controls   = counts['up_controls']
    down_controls = counts['down_controls']
    n_cases       = counts['n_cases']
    n_controls    = counts['n_controls']

    # NRI components (Pencina et al. 2008, Stat Med 27:157-172)
    nri_events     = (up_cases      - down_cases)    / n_cases
    nri_non_events = (down_controls - up_controls)   / n_controls
    nri            = nri_events + nri_non_events

    # Standard error: SE(NRI+) = sqrt((up + down) / n^2) = sqrt(up + down) / n
    se_events     = np.sqrt(up_cases    + down_cases)    / n_cases
    se_non_events = np.sqrt(up_controls + down_controls) / n_controls
    se_nri        = np.sqrt(se_events**2 + se_non_events**2)

    z_stat  = nri / se_nri if se_nri > 0 else np.nan
    p_value = 2 * (1 - norm.cdf(abs(z_stat))) if not np.isnan(z_stat) else np.nan

    ci_low  = nri - 1.96 * se_nri
    ci_high = nri + 1.96 * se_nri

    print(f"NRI (Events): {nri_events:.3f}")
    print(f"  Cases: {up_cases} moved up, {down_cases} moved down")
    print(f"NRI (Non-events): {nri_non_events:.3f}")
    print(f"  Controls: {down_controls} moved down, {up_controls} moved up")
    print(f"\nTotal NRI: {nri:.3f}")
    print(f"95% CI: [{ci_low:.3f}, {ci_high:.3f}]")
    print(f"z-statistic: {z_stat:.3f}")
    print(f"P-value: {p_value:.4e}")

    if p_value < 0.05:
        if nri > 0:
            print(f"   -> {prs2_name} provides significant improvement in reclassification")
        else:
            print(f"   -> {prs1_name} provides better reclassification")
    else:
        print(f"   -> No significant difference in reclassification")

    # IDI: improvement in discrimination slope (Pencina 2008)
    # IDI = (mean_cases_new - mean_controls_new) - (mean_cases_old - mean_controls_old)
    mean_prs1_cases    = data.loc[is_case,  f'scaled_prs_{prs1_name}'].mean()
    mean_prs1_controls = data.loc[~is_case, f'scaled_prs_{prs1_name}'].mean()
    mean_prs2_cases    = data.loc[is_case,  f'scaled_prs_{prs2_name}'].mean()
    mean_prs2_controls = data.loc[~is_case, f'scaled_prs_{prs2_name}'].mean()

    idi = (mean_prs2_cases - mean_prs2_controls) - (mean_prs1_cases - mean_prs1_controls)
    print(f"\nIDI (Integrated Discrimination Improvement): {idi:.4f}")

    return {
        'test':           'NRI',
        'statistic':      nri,
        'pvalue':         p_value,
        'conf_int':       (ci_low, ci_high),
        'nri_events':     nri_events,
        'nri_non_events': nri_non_events,
        'idi':            idi,
    }


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
    df["high_prs1"] = df[prs_method1]
    
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
            
            if clinical_measure != 'None':
                #get the people with low clinical measures
                df = dfFull[dfFull[measure] == 9]

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
        
    
    # Calculate percentage of high risk
    high_risk_percent = (df[binary_col] == 1).mean() * 100
    
    plt.title(f'Distribution of {measure}\nThreshold: {threshold:.2f} ({percentile:.1f}th percentile, {high_risk_percent:.1f}% classified as high risk)')
    plt.xlabel(measure)
    plt.ylabel('Frequency')
    plt.grid(alpha=0.3)
    
    return plt.gcf()



def compare_prs_performance(df, clinical_measures, figPath, file_ext, prs_methods, risk_thresholds=None,
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
    prs_methods : list of columns to inlcude using different methods
        scaled_prs or centile_bin
        
    Returns:
    --------
    dict
        Dictionary containing performance metrics
    pandas.DataFrame
        DataFrame with added binary columns
    """
    results = {}

    
    df_copy = df.copy()
    # Create list of binary clinical measures
    binary_measures = [f"{measure}_binary" for measure in clinical_measures 
                       if f"{measure}_binary" in df_copy.columns]
    
    # If no binary measures were created, use original measures (assuming they're already binary)
    if not binary_measures:
        binary_measures = clinical_measures
        

    
    resultsLowClinicalHighPRS = pd.DataFrame()
    resultsLowClinical = pd.DataFrame()

    
    
    # Calculate calibration
#   results_calibration = calculate_calibration(df_binary, prs_methods, binary_measures,'PHENOTYPE')
#   results['calibration'] = results_calibration
    #results['auc_plots'] = auc_plots
    
    
    results['nri'] = {}
    for i, prs1 in enumerate(prs_methods):
        if 'scaled_prs' in prs1:
            threshold = df_copy[prs1].quantile(0.80)
            df_copy[f'{prs1}_binary'] = (df_copy[prs1] >= threshold).astype(int)
            prs_col1 = f'{prs1}_binary'
        else:
            if prs1.startswith('bin_'):
                base_prs1 = prs1[4:]                       # 'bin_main' → 'main'
            else:
                base_prs1 = prs1.replace('_centile_bin', '')  # old naming, backward compat
            prs_col1 = f'{base_prs1}_high_risk'
        
        results['nri'][prs1] = {}
        for measure in binary_measures:
            results['nri'][f"{prs1}_vs_{measure}"] = {}
            try:
                nri, nri_events, nri_non_events = calculate_nri(
                    df_copy, measure, prs_col1, True
                )
                
                
                results['nri'][f"{prs1}_vs_{measure}"][measure] = {
                    'nri': nri,
                    'nri_events': nri_events,
                    'nri_non_events': nri_non_events
                }

            except Exception as e:
                results['nri'][f"{prs1}_vs_{measure}"][measure] = {'error': str(e)}
            
            try:
                #calculate nri with PRS_combined and low clinical measure
#               results['nri'][f"{prs1}_vs_{measure}"] = {}
                nri, nri_events, nri_non_events = calculate_nri(
                    df_copy[df_copy[measure] == 0], measure, prs_col1, True
                )
                
                
                results['nri'][f"{prs1}_vs_{measure}"][f"{measure}_low"] = {
                    'nri': nri,
                    'nri_events': nri_events,
                    'nri_non_events': nri_non_events
                }
            
            except Exception as e:
                results['nri'][f"{prs1}_vs_{measure}"][f"{measure}_low"] = {'error': str(e)}
                
                
            fig = plot_reclassification_table(df_copy, measure, prs_col1)
            fig.savefig(f'{figPath}/reclassificationHeatMap.{measure}v{prs1}{file_ext}.png',dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            fig = plot_low_clinical_reclassification_table(df_copy, measure, prs_col1)
            fig.savefig(f'{figPath}/reclassificationHeatMap.LowOnly_{measure}v{prs1}{file_ext}.png',dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            
        
    #       if prs1 == 'combined_centile_bin':
            fig,nri_summary = plot_nri_from_reclassification(df_copy,binary_measures,prs_col1)
            plt.savefig(f'{figPath}/nri_comparison_clinical_vars_bar_chart.{prs1}.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            
        # Calculate AUC for all methods and binary measures
        resultsLowClinicalTemp, resultsLowClinicalHighPRSTemp = calculate_auc(df_copy, prs1, prs_col1, binary_measures)
        
        resultsLowClinical = pd.concat([resultsLowClinical,resultsLowClinicalTemp],ignore_index=False)
        resultsLowClinicalHighPRS = pd.concat([resultsLowClinicalHighPRS,resultsLowClinicalHighPRSTemp],ignore_index=False)
        
        
        results['auc_low_clinical'] = resultsLowClinical
        results['auc_low_clinical_high_prs'] = resultsLowClinicalHighPRS 
        
        for j, prs2 in enumerate(prs_methods):
            if i >= j:  # Skip self-comparisons and repeated comparisons
                continue
            try:
                if 'scaled_prs' in prs2:
                    threshold = df_copy[prs2].quantile(0.80)
                    df_copy[f'{prs2}_binary'] = (df_copy[prs2] >= threshold).astype(int)
                    prs_col1 = f'{prs2}_binary'

                else:
                    base_prs2 = prs2.replace('_centile_bin','')
                    prs_col2 = f'{base_prs2}_high_risk'
        
                key = f"{prs1}_vs_{prs2}"
                results['nri'][key] = {}
    

                nri, nri_events, nri_non_events = calculate_nri(
                    df_copy, prs_col1, prs_col2, False
                )
                results['nri'][key][f'{measure}_low'] = {
                    'nri': nri,
                    'nri_events': nri_events,
                    'nri_non_events': nri_non_events
                }
                
                fig = plot_reclassification_table(df_copy, prs_col1, prs_col2)
                fig.savefig(f'{figPath}/reclassificationHeatMap.{prs1}v{prs2}{file_ext}.png',dpi=150, bbox_inches='tight')
                plt.close(fig)

            except Exception as e:
                results['nri'][key][f'{measure}_low'] = {'error': str(e)}
    
    return results



def impute_clinical_data(train_df, test_df, clinical_columns, 
                        outcome_column,
                        method='mean', visualize=True):
    """
    Impute missing values in clinical measures for PRS analysis
    
    Parameters:
    -----------
    train_df : pandas.DataFrame
        DataFrame containing clinical measures,  and outcome variables for training data
    test_df : pandas.DataFrame
        DataFrame containing clinical measures, and outcome variables for training data
    clinical_columns : list
        List of clinical measure column names that might contain missing values
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
    predictor_cols = clinical_columns.copy()
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

    
    
def run_clinical_analysis(pheno, pheno_data, results_path, clinical_source_path=None):
    """
    Run clinical-measure performance analysis on the combined holdout PRS set.

    Parameters
    ----------
    pheno : str
        Phenotype name (e.g. 'type2Diabetes').
    pheno_data : str
        Root of the combinedAnalysis directory (contains scores/, figures/).
    results_path : str
        Broader results directory containing participant_environment.csv.
    clinical_source_path : str, optional
        Directory that contains clinicalEnvironmentHoldout.csv.  When None
        (default) the function searches for the file in the sibling
        productEpi/ and summedEpi/ sub-directories of the parent of
        pheno_data, since combinedAnalysis does not have its own copy.
    """
    # Generate simulated data
    figPath = f'{pheno_data}/figures/clinicalFigures'
    scoresPath = f'{pheno_data}/scores/clinicalMeasures'
    # Creates directory and all parent directories if needed
    # Does nothing if it already exists (no error)
    os.makedirs(f'{figPath}', exist_ok=True)
    
    os.makedirs(f'{scoresPath}', exist_ok=True)
    
    df = pd.read_csv(f'{pheno_data}/scores/combinedPRSGroups.filtered.csv')
    
    #get the filtered models to use in analysis
#   prs_columns = [m for m in df.columns if 'bin' in m]

    validation_df = pd.read_csv(f'{pheno_data}/scores/combinedPRSGroups.holdout.filtered.csv')
    
    #prs columns to use in analysis
#   prs_columns = [col for col in df.columns if 'scaled_prs' in col]
    #prs_columns = [col for col in prs_columns if 'scaled' not in col]
    
    if pheno == 'type2Diabetes':
#       clinical_columns = ['Glycated haemoglobin (HbA1c)','Body mass index (BMI)','Glucose']
        # Specify custom thresholds for some measures (optional)
        clinical_thresholds = {
            'Glycated haemoglobin (HbA1c)': 41,
            'Body mass index (BMI)': 25,
            'Glucose' : 12
        }
        
    else:
#       clinical_columns = ['Basal metabolic rate','Urea','Haemoglobin concentration']
        clinical_thresholds = None
    
    # Resolve the clinical column source file.
    # combinedAnalysis does not contain clinicalEnvironmentHoldout.csv;
    # it is identical across productEpi and summedEpi, so we locate it
    # in one of those sibling directories when not supplied explicitly.
    if clinical_source_path is None:
        parent = os.path.dirname(os.path.abspath(pheno_data))
        for candidate in ('productEpi', 'summedEpi'):
            candidate_file = os.path.join(parent, candidate, 'clinicalEnvironmentHoldout.csv')
            if os.path.exists(candidate_file):
                clinical_source_path = candidate_file
                print(f'  [clinical] Using clinicalEnvironmentHoldout.csv from {candidate}/')
                break
        if clinical_source_path is None:
            raise FileNotFoundError(
                'clinicalEnvironmentHoldout.csv not found in productEpi/ or summedEpi/ '
                f'siblings of {pheno_data}.  Pass clinical_source_path= explicitly.'
            )
    clinical_columns = pd.read_csv(clinical_source_path, nrows=1)
    clinical_columns = clinical_columns.set_index('IID').columns.tolist()
        
    clinical_data = pd.read_csv(f'{results_path}/participant_environment.csv',usecols=['Participant ID']+clinical_columns)
    clinical_data.rename(columns={'Participant ID':'IID'},inplace=True)
    
    #combine clinical data with PHENOTYPE for imputation methods
    data = df[['IID','PHENOTYPE']].merge(clinical_data, on=['IID'],how='left')
    validation_data = validation_df[['IID','PHENOTYPE']].merge(clinical_data, on=['IID'],how='left')
    
    
    train_imputed, test_imputed = impute_clinical_data(data, validation_data, clinical_columns, 'PHENOTYPE', method='mean', visualize=False)
    
    #merge with clinical train_imputed 
    train_data_clinical = df.merge(train_imputed, on=['IID','PHENOTYPE'],how='left')
    
    # Update values
    train_data_clinical.loc[train_data_clinical["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
    train_data_clinical.loc[train_data_clinical["PHENOTYPE"] == 2, "PHENOTYPE"] = 1
    
    holdout_data_clinical = validation_df.merge(test_imputed, on=['IID','PHENOTYPE'],how='left')
    holdout_data_clinical.loc[holdout_data_clinical["PHENOTYPE"] == 1, "PHENOTYPE"] = 0
    holdout_data_clinical.loc[holdout_data_clinical["PHENOTYPE"] == 2, "PHENOTYPE"] = 1

    
    # First convert clinical measures to binary if needed
    df_binary_test, df_binary_validation, used_thresholds, determined_directions = convert_to_binary(
        train_data_clinical, holdout_data_clinical, clinical_columns, thresholds=clinical_thresholds, 
        high_risk_quintile=True, risk_directions=None,outcome_column='PHENOTYPE'
    )
    
    
#   for item_tuple in [(df_binary_validation,'holdout'),(df_binary_test,'')]:
#   for item_tuple in [(df_binary_validation,'holdout')]:
            
        
    # Print thresholds and risk directions used
    df_binary = df_binary_validation.copy()
    file_ext = 'holdout'
    
    #create_sankey_plot_clinical_data(df_binary,figPath,use_epi_main=False)
    
    for eval_type in ['bin', 'scaled_prs']:
        if eval_type == 'bin':
            # Use startswith — 'bin' substring would also match _binary clinical columns
            prs_methods = [
                col for col in df_binary
                if col.startswith('bin_') or col.endswith('_centile_bin')
            ]
        else:
            prs_methods = [
                col for col in df_binary
                if 'scaled_prs' in col and 'threshold_value' not in col and 'any' not in col
            ]
    
        # Run comprehensive comparison with direction-aware threshold conversion
        results = compare_prs_performance(
            df_binary, clinical_columns, figPath, file_ext, prs_methods,
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
        auc_low_clinical = results['auc_low_clinical']
#       results['auc_low_clinical'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossPRSLowClinicalMeasure{file_ext}.csv')
        
        # Print AUC results
        print("\nAUC Results low clinical high prs")
        print(results['auc_low_clinical_high_prs'])
        auc_low_clinical_high_prs = results['auc_low_clinical_high_prs']
#       results['auc_low_clinical_high_prs'].to_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/wrtClinical/aucAcrossHighPRSLowClinicalMeasure{file_ext}.csv')
        
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
        
        #clinical_auc.to_csv(f"{scoresPath}/prsAUC.wrt.clinicalMeasures{file_ext}.csv",index=False)
        

        
        #put nri data into a csv
        df_nri = pd.DataFrame()
        for prs in results['nri'].keys():
            prs_nri = pd.DataFrame(results['nri'][prs]).T
            prs_nri['prs_calc'] = prs
            prs_nri.reset_index(inplace=True)
            prs_nri.columns = ['clinical_measure'] + list(prs_nri.columns[1:])
            df_nri = pd.concat([df_nri,prs_nri],ignore_index=True)
            #df_nri.to_csv(f"{scoresPath}/nri.wrt.clinicalMeasures{file_ext}.csv",index=False)

        
        
        for k in ['risk_directions','thresholds_used']:
            try:
                thresholds_risk[k] = pd.DataFrame(results[k],index=[k]).T
            except UnboundLocalError:
                thresholds_risk = pd.DataFrame(results[k],index=[k]).T
                
        thresholds_risk.reset_index(inplace=True)
        thresholds_risk.columns = ['clinical_measure'] + list(thresholds_risk.columns[1:])
        
        #thresholds_risk.to_csv(f"{scoresPath}/thresholdsRiskDirection.wrt.clinicalMeasures{file_ext}.csv",index=False)
        
    # Export results to Excel with multiple sheets
        with pd.ExcelWriter(f'{scoresPath}/clinicalMeasurePerformanceResults{file_ext}.{eval_type}.xlsx') as writer:
            thresholds_risk.to_excel(writer, sheet_name='thresholds RiskDirection')
            df_nri.to_excel(writer, sheet_name='nri')
            clinical_auc.to_excel(writer, sheet_name='auc ')
            auc_low_clinical_high_prs.to_excel(writer, sheet_name='auc low clinical high prs')
            auc_low_clinical.to_excel(writer, sheet_name='auc low clinical')
            
        print(f"\nResults exported to 'clinicalMeasurePerformanceResults{file_ext}.{eval_type}.xlsx'")
        
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
        fig.savefig(f'{figPath}/riskDistribution.{measure}{file_ext}.png')
        plt.close(fig)
        
    
    #save binary file with prs calculation and clinical marker data
    df_binary.to_csv(f"{scoresPath}/combinedPRS.holdout.ClinicalMeasures.csv",index=False)
        
        
        
        
        

# Backward-compatible alias so existing callers using `main` still work
main = run_clinical_analysis


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Calculating clinical performance measures against prs ....")
    parser.add_argument("--pheno_data",help="path to pheno results directory")
    parser.add_argument("--pheno",help="path to pheno")
    parser.add_argument("--results_path",help="path to cleaned covariate file")
    
    args = parser.parse_args()
    
    pheno = args.pheno or os.environ.get("PHENO")
    print(f"[PYTHON] Phenotype is : {pheno}")
    
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
    results_path = args.results_path or os.environ.get("RESULTS_PATH")
    print(f"[PYTHON] Reading from: {results_path}")
    
#   pheno = "type2Diabetes"    
#   pheno_data = "/Users/kerimulterer/prsInteractive/results/type2Diabetes/combinedAnalysis"
#   results_path = "/Users/kerimulterer/prsInteractive/results"
    
    
    if not pheno_data:
        raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
    
    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
        
    if not results_path:
        raise ValueError("You must provide a results path via --results_path or set the RESULTS_PATH environment variable.")
        
    run_clinical_analysis(pheno, pheno_data, results_path)