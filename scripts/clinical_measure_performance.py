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

def create_sankey_plot_clinical_data(df,figPath,use_epi_main=False):
    
    if use_epi_main:
        prs_mathods = ['main','epi+main','epi','cardio']
    else:
        prs_methods = ['main', 'epi', 'cardio'] 
        
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
                node_colors.append(COHORT_COLORS['cardio'])
                node_y.append(0.1)  # Top
            elif 'EPI' in label and 'main' not in label.lower():
                node_colors.append(COHORT_COLORS['epi'])
                node_y.append(0.4)  # Middle
            elif 'MAIN' in label:
                node_colors.append(COHORT_COLORS['main'])
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
    
    Clinical measures: 1=Low Risk, 0=High Risk
    PRS: 1=High Risk, 0=Low Risk
    
    When prs_col contains 'combined', creates stacked bars showing the contribution
    of each cohort (main, epi, cardio) to NRI for cases and controls separately.
    
    For each clinical variable:
    - Controls: Count those moved from high risk (0) to low risk (0) by PRS
        (i.e., clinical says high risk [0], PRS correctly says low risk [0])
    - Cases: Count those moved from low risk (1) to high risk (1) by PRS
        (i.e., clinical says low risk [1], PRS correctly says high risk [1])
    
    Parameters:
    -----------
    df : pandas.DataFrame
            DataFrame with 'PHENOTYPE' column (1=case, 0=control)
    clinical_vars : list of clinical var str
            list clinical risk category column names (1=low risk, 0=high risk)
    prs_col : str
            Column name for PRS risk categories (1=high risk, 0=low risk)
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
            df.loc[controls, risk_cat1],  # Clinical risk (rows): 1=low, 0=high
            df.loc[controls, prs_col],     # PRS risk (columns): 1=high, 0=low
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
        total_controls_high = reclass_controls.loc[0].sum() if 1 in reclass_controls.index else 0
        
        # NRI- for controls (proportion correctly kept at/moved to low risk)
        nri_controls = controls_correct / len(df[controls]) if len(df[controls]) > 0 else 0
        
        # =====================================================================
        # CASES: Clinical says LOW risk (1), want PRS to correctly say HIGH risk (1)
        # =====================================================================
        reclass_cases = pd.crosstab(
            df.loc[cases, risk_cat1],      # Clinical risk (rows): 1=low, 0=high
            df.loc[cases, prs_col],         # PRS risk (columns): 1=high, 0=low
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
        total_cases_low = reclass_cases.loc[1].sum() if 0 in reclass_cases.index else 0
        
        # NRI+ for cases (proportion correctly moved to high risk)
        nri_cases = cases_correct / len(df[cases]) if len(df[cases]) > 0 else 0
        
        # =====================================================================
        # COHORT-SPECIFIC NRI (only for combined PRS)
        # =====================================================================
        cohort_nri = {}
        if is_combined:
            
            
            # Derive individual PRS column names
            main_col = prs_col.replace('combined', 'main')
            epi_col = prs_col.replace('combined', 'epi')
            cardio_col = prs_col.replace('combined', 'cardio')
            
            # For CASES: Clinical LOW (0) → PRS HIGH (1)
            # Identify which cohort flagged each correctly reclassified case
            cases_reclass_mask = (df[risk_cat1] == 0) & (df[prs_col] == 1) & cases
            
            # Priority: main > epi > cardio
            cases_main = ((df[main_col] == 1) & cases_reclass_mask).sum()
            cases_epi = ((df[epi_col] == 1) & cases_reclass_mask & (df[main_col] == 0)).sum()
            cases_cardio = ((df[cardio_col] == 1) & cases_reclass_mask & 
                           (df[main_col] == 0) & (df[epi_col] == 0)).sum()
            
            # For CONTROLS: Clinical HIGH (1) → PRS LOW (0)
            # These controls are NOT flagged by any PRS (all PRSs agree = low risk)
            # We'll distribute them proportionally based on overall cohort sizes
            controls_reclass_mask = (df[risk_cat1] == 1) & (df[prs_col] == 0) & controls
            controls_reclass_total = controls_reclass_mask.sum()
            
            # Get overall cohort sizes for proportional attribution
            total_main = (df[main_col] == 1).sum()
            total_epi = (df[epi_col] == 1).sum()
            total_cardio = (df[cardio_col] == 1).sum()
            total_all = total_main + total_epi + total_cardio
            
            if total_all > 0:
                controls_main = controls_reclass_total * (total_main / total_all)
                controls_epi = controls_reclass_total * (total_epi / total_all)
                controls_cardio = controls_reclass_total * (total_cardio / total_all)
            else:
                controls_main = controls_epi = controls_cardio = 0
                
            # Calculate NRI components by cohort
            n_cases = len(df[cases])
            n_controls = len(df[controls])
            
            cohort_nri = {
                'cases_main': cases_main / n_cases if n_cases > 0 else 0,
                'cases_epi': cases_epi / n_cases if n_cases > 0 else 0,
                'cases_cardio': cases_cardio / n_cases if n_cases > 0 else 0,
                'cases_main_n': cases_main,
                'cases_epi_n': cases_epi,
                'cases_cardio_n': cases_cardio,
                'controls_main': controls_main / n_controls if n_controls > 0 else 0,
                'controls_epi': controls_epi / n_controls if n_controls > 0 else 0,
                'controls_cardio': controls_cardio / n_controls if n_controls > 0 else 0,
                'controls_main_n': int(controls_main),
                'controls_epi_n': int(controls_epi),
                'controls_cardio_n': int(controls_cardio)
            }
            
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
        # STACKED BARS: Cases and Controls NRI by cohort
        # =====================================================================
        width = 0.35
        
#       cohort_colors = {
#           'main': '#E74C3C',      # Red
#           'epi': '#3498DB',       # Blue
#           'cardio': '#F39C12'     # Orange
#       }
        
        # CASES STACKED BAR (left)
        bars_cases_main = ax.bar(
            x - width/2,
            nri_summary['cases_main'],
            width,
            label=f'Main PRS',
            color=COHORT_COLORS['main'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
        bars_cases_epi = ax.bar(
            x - width/2,
            nri_summary['cases_epi'],
            width,
            bottom=nri_summary['cases_main'],
            label=f'Epistatic PRS',
            color=COHORT_COLORS['epi'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
        bars_cases_cardio = ax.bar(
            x - width/2,
            nri_summary['cases_cardio'],
            width,
            bottom=nri_summary['cases_main'] + nri_summary['cases_epi'],
            label=f'Cardio PRS',
            color=COHORT_COLORS['cardio'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
        # CONTROLS STACKED BAR (right)
        bars_controls_main = ax.bar(
            x + width/2,
            nri_summary['controls_main'],
            width,
            color=COHORT_COLORS['main'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
        bars_controls_epi = ax.bar(
            x + width/2,
            nri_summary['controls_epi'],
            width,
            bottom=nri_summary['controls_main'],
            color=COHORT_COLORS['epi'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
        
        bars_controls_cardio = ax.bar(
            x + width/2,
            nri_summary['controls_cardio'],
            width,
            bottom=nri_summary['controls_main'] + nri_summary['controls_epi'],
            color=COHORT_COLORS['cardio'],
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5
        )
        
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
            
        # Add "Cases" and "Controls" labels
        y_pos = ax.get_ylim()[1] * 0.95
#       ax.text(x[len(x)//2] - width/2, y_pos, 'Cases (Clinical Low→PRS High)',
#              ha='center', va='top', fontsize=11, fontweight='bold',
#              bbox=dict(boxstyle='round', facecolor=colors['cases'], alpha=0.3))
#       ax.text(x[len(x)//2] + width/2, y_pos, 'Controls (Clinical High→PRS Low)',
#              ha='center', va='top', fontsize=11, fontweight='bold',
#              bbox=dict(boxstyle='round', facecolor=colors['controls'], alpha=0.3))
        
    else:
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


def calculate_auc(dfFull, prs_method, clinical_measures):
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
            
        try:
            if 'scaled_prs' in prs_method:
                threshold = df[prs_method].quantile(.80)
                df['high_prs'] = (df[prs_method] >= threshold).astype(int)
            else:
                df['high_prs'] = (df[prs_method] >= 800).astype(int)
            dfHigh = df[df['high_prs'] == 1]
            auc = roc_auc_score(dfHigh['PHENOTYPE'], dfHigh[prs_method])
            resultsTop20PRS.loc[prs_method, measure] = auc
            
        except Exception as e:
            resultsTop20PRS.loc[prs_method, measure] = np.nan
            print(f"Error calculating AUC for high risk {prs_method} vs low risk {measure}: {e}")
                
    return resultsAllPRS, resultsTop20PRS


def calculate_nri(df, prs_method1, prs_method2,clinical=True):
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
    # Categorize: Top 20% as "High Risk", others as "Low Risk"
    # if using centile_bin, high risk column already created

    # Extract events and non-events
    events = df['PHENOTYPE'] == 1
    non_events = df['PHENOTYPE'] == 0
    
    if clinical:
        # Ensure the clinical measure is binary
        if not set(df[prs_method1].unique()).issubset({0, 1}):
            raise ValueError(f"Clinical measure {clinical_measure} must be binary (0/1)")

#           
#       # Extract events and non-events
#       events = df[prs_method1] == 1
#       non_events = df[prs_method1] == 0
#       
#   else:
#       # Calculate NRI for all pairs of PRS methods
        
    risk_cat1 = df[prs_method1]
    risk_cat2 = df[prs_method2]
    

    
    # Calculate proportions of up-classifications and down-classifications
    
    #cases going from low risk clinical to high risk prs (low risk clinical = 0, high risk prs = 1)
    up_events = np.mean(risk_cat2[events] > risk_cat1[events])
    #cases going from high_risk prs to low-risk clinical measure
    down_events = np.mean(risk_cat2[events] < risk_cat1[events])
    
    #controls going from high risk clinical to low risk prs (low risk clinical = 0, high risk prs = 1)
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
            df_copy[f'{prs1}_binary'] = df[prs1].quantile(0.8).astype(int)
            prs_col1 = f'{prs1}_binary'
        else:
            base_prs1 = prs1.replace('_centile_bin','')
            prs_col1 = f'{base_prs1}_high_risk'
        
        results['nri'][prs1] = {}
        for measure in binary_measures:
            try:
                results['nri'][f"{prs1}_vs_{measure}"] = {}
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
        resultsLowClinicalTemp, resultsLowClinicalHighPRSTemp = calculate_auc(df_copy, prs1, binary_measures)
        
        resultsLowClinical = pd.concat([resultsLowClinical,resultsLowClinicalTemp],ignore_index=False)
        resultsLowClinicalHighPRS = pd.concat([resultsLowClinicalHighPRS,resultsLowClinicalHighPRSTemp],ignore_index=False)
        
        for j, prs2 in enumerate(prs_methods):

            if i >= j:  # Skip self-comparisons and repeated comparisons
                continue
            
            if 'scaled_prs' in prs2:
                df_copy[f'{prs2}_binary'] = df[prs2].quantile(0.8).astype(int)
                prs_col2 = f'{prs2}_binary'
            else:
                base_prs2 = prs2.replace('_centile_bin','')
                prs_col2 = f'{base_prs2}_high_risk'
    
            key = f"{prs1}_vs_{prs2}"
            results['nri'][key] = {}
    
    
            #get the calibration key in for comparing prs calculations
#                   results['calibration'][key] = {}
    
            try:
                if 'scaled_prs' in prs2:
                    df_copy[f'{prs1}_binary'] = df[prs1].quantile(0.8).astype(int)
                    df_copy[f'{prs2}_binary'] = df[prs2].quantile(0.8).astype(int)
                    prs_col1 = f'{prs1}_binary'
                    prs_col2 = f'{prs2}_binary'
                    
                else:
                    base_prs1 = prs1.replace('_centile_bin','')
                    prs_col1 = f'{base_prs1}_high_risk'
                    base_prs2 = prs2.replace('_centile_bin','')
                    prs_col2 = f'{base_prs2}_high_risk'

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
        
    results['auc_low_clinical'] = resultsLowClinical
    results['auc_low_clinical_high_prs'] = resultsLowClinicalHighPRS    
    
    
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

    
    
def main(pheno,pheno_data,results_path):
    """Run an example using simulated data to demonstrate the function"""
    # Generate simulated data
    figPath = f'{pheno_data}/figures/clinicalFigures'
    scoresPath = f'{pheno_data}/scores/clinicalMeasures'
    # Creates directory and all parent directories if needed
    # Does nothing if it already exists (no error)
    os.makedirs(f'{figPath}', exist_ok=True)
    
    os.makedirs(f'{scoresPath}', exist_ok=True)
    
    df = pd.read_csv(f'{pheno_data}/scores/combinedPRSGroups.csv')
    validation_df = pd.read_csv(f'{pheno_data}/scores/combinedPRSGroups.holdout.csv')
    
    #prs columns to use in analysis
#   prs_columns = [col for col in df.columns if 'scaled_prs' in col]
    #prs_columns = [col for col in prs_columns if 'scaled' not in col]
    
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
    for item_tuple in [(df_binary_validation,'holdout')]:
            
        
    # Print thresholds and risk directions used
        df_binary = item_tuple[0].copy()
        file_ext = item_tuple[1]
        
        create_sankey_plot_clinical_data(df_binary,figPath,use_epi_main=False)
        
        for eval_type in ['centile_bin','scaled_prs']:   
            prs_methods = [col for col in df_binary if eval_type in col]
        
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
            
        #save binary file with prs calculation and clinical marker data
        df_binary.to_csv(f"{scoresPath}/combinedPRS.holdout.ClinicalMeasures.csv",index=False)
        
        
        
        
        
        
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
    
    pheno = "type2Diabetes"    
    pheno_data = "/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi"
    results_path = "/Users/kerimulterer/prsInteractive/results"
    
    
    if not pheno_data:
        raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
    
    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
        
    if not results_path:
        raise ValueError("You must provide a results path via --results_path or set the RESULTS_PATH environment variable.")
        
    main(pheno,pheno_data,results_path)
    

    
    
    

    

    
    
    