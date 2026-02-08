#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
from helper.draw_plots import *
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint

labels = {'main': 'PRS$_G$', 'epi': 'PRS$_{G×G}$', 
                    'cardio': 'PRS$_{G×G×E}$', 'combined': 'PRS$_{comb}$'}
                
def plot_odds_ratios_comparison(df, percentiles=[1,3,5,10,20],prs_methods=None, figsize=(12, 8),
                                    show_significance=True, reference_line=True,
                                    save_path=None, dpi=300):
        """
        Create scatter plot showing odds ratios with error bars for different PRS methods
        across percentile thresholds.
        
        Parameters:
        -----------
        df : DataFrame
                Data with columns: percentiles, odds_ratios, p_values, CI_Upper, CI_lower, model
        prs_methods : list
                List of PRS methods to include (if None, uses all in df)
        figsize : tuple
                Figure dimensions
        show_significance : bool
                Whether to mark significant results
        reference_line : bool
                Whether to show OR=1 reference line
        save_path : str
                Path to save figure
        dpi : int
                Resolution
        
        Returns:
        --------
        fig, ax : matplotlib figure and axis objects
        """
    
        # Filter to specified PRS methods if provided
        if prs_methods is not None:
                df = df[df['model'].isin(prs_methods)].copy()
                df = df[~df['model'].str.contains('centile')]
            
        # Get unique percentiles and models
        models = df['model'].unique()
        
        if percentiles is not None:
            df = df[df['percentiles'].isin(percentiles)]
        else:
            percentiles = df['percentiles'].unique()
    
        # Sort percentiles if they follow a pattern (e.g., "Top 5%", "Top 10%")
        try:
            percentiles = sorted(percentiles, key=lambda x: float(x.replace('Top ', '').replace('%', '')))
        except:
            pass
            
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
    
        # Calculate positions for each model within each percentile
        n_models = len(models)
        group_width = 0.8
        bar_width = group_width / n_models
    
        # Create x positions for percentiles
        x_positions = np.arange(len(percentiles))
    
        # Plot each model
        for model_idx, model in enumerate(models):
                # Filter data for this model
                model_data = df[df['model'] == model].copy()
            
                # Sort by percentile order
                model_data['percentile_order'] = model_data['percentiles'].map(
                        {p: i for i, p in enumerate(percentiles)}
                )
                model_data = model_data.sort_values('percentile_order')
                model_data.dropna(inplace=True)
            
                # Calculate offset for this model
                offset = (model_idx - n_models/2 + 0.5) * bar_width
                x_pos = x_positions + offset
            
                # Get colors and labels
                color = COHORT_COLORS.get(model, 'gray')
                label = labels.get(model, model)
            
                # Calculate error bars
                yerr_lower = model_data['odds_ratios'] - model_data['CI_lower']
                yerr_upper = model_data['CI_Upper'] - model_data['odds_ratios']
            
                # Plot scatter with error bars
                ax.errorbar(x_pos, model_data['odds_ratios'], 
                                        yerr=[yerr_lower, yerr_upper],
                                        fmt='o', markersize=8, capsize=5, capthick=2,
                                        color=COHORT_COLORS[model], label=labels[model], linewidth=2, alpha=0.8)
            
                # Mark significant points if requested
#               if show_significance:
#                       sig_mask = model_data['p_values'] < 0.05
#                       if sig_mask.any():
#                               ax.scatter(x_pos[sig_mask], model_data.loc[sig_mask, 'odds_ratios'],
#                                                   marker='*', s=200, color=COHORT_COLORS[model], edgecolors='black',
#                                                   linewidths=1, zorder=5)
                            
        # Add reference line at OR = 1
        if reference_line:
                ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, 
                                    alpha=0.7, label='OR = 1 (No effect)', zorder=1)
            
        # Formatting
        ax.set_xlabel('Risk Percentile Threshold', fontsize=13, fontweight='bold')
        ax.set_ylabel('Odds Ratio (OR)', fontsize=13, fontweight='bold')
        ax.set_title('T2D Risk: Odds Ratios by PRS Model and Percentile Threshold', 
                                fontsize=14, fontweight='bold', pad=20)
    
        # Set x-ticks
        ax.set_xticks(x_positions)
        ax.set_xticklabels(percentiles, fontsize=11)
    
        # Y-axis formatting
        ax.set_ylim(bottom=0)
        ax.tick_params(axis='y', labelsize=11)
    
        # Legend
#       legend = ax.legend(loc='upper left', frameon=True, fancybox=True, 
#                                           shadow=True, fontsize=11)
    
#       if show_significance:
#               # Add note about significance markers
#               ax.text(0.98, 0.02, '* p < 0.05', transform=ax.transAxes,
#                               fontsize=10, ha='right', va='bottom',
#                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
            
        # Grid
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
    
        plt.tight_layout()
    
        if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
                print(f"Figure saved to {save_path}")
            
        return fig, ax


def plot_odds_ratios_line_plot(df, prs_methods=None, figsize=(12, 8),
                                                                show_ci_bands=True, save_path=None, dpi=300):
        """
        Create line plot showing odds ratios with confidence interval bands across percentiles.
        
        Parameters:
        -----------
        df : DataFrame
                Data with columns: percentiles, odds_ratios, p_values, CI_Upper, CI_lower, model
        prs_methods : list
                List of PRS methods to include
        figsize : tuple
                Figure dimensions
        show_ci_bands : bool
                Whether to show CI as shaded bands
        save_path : str
                Path to save figure
        dpi : int
                Resolution
        
        Returns:
        --------
        fig, ax : matplotlib figure and axis objects
        """
    
        # Filter to specified PRS methods if provided
        if prs_methods is not None:
                df = df[df['model'].isin(prs_methods)].copy()
            
    
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
    
        # Get unique models
        models = df['model'].unique()
    
        # Extract numeric percentile values for x-axis
        df['percentile_num'] = df['percentiles'].str.extract(r'(\d+\.?\d*)').astype(float)
    
        # Plot each model
        for model in models:
                model_data = df[df['model'] == model].sort_values('percentile_num')
            
                color = COHORT_COLORS.get(model, 'gray')
                label = labels.get(model, model)
            
                # Plot line
                ax.plot(model_data['percentile_num'], model_data['odds_ratios'],
                                marker='o', markersize=8, linewidth=2.5, color=color, label=label)
            
                # Add CI bands if requested
                if show_ci_bands:
                        ax.fill_between(model_data['percentile_num'],
                                                        model_data['CI_lower'],
                                                        model_data['CI_Upper'],
                                                        color=color, alpha=0.2)
                    
        # Reference line
        ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, 
                            alpha=0.7, label='OR = 1 (No effect)')
    
        # Formatting
        ax.set_xlabel('Risk Percentile Threshold (%)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Odds Ratio (OR)', fontsize=13, fontweight='bold')
        ax.set_title('T2D Risk: Odds Ratios by PRS Model Across Percentile Thresholds', 
                                fontsize=14, fontweight='bold', pad=20)
    
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=11)
        ax.grid(True, alpha=0.3)
    
        plt.tight_layout()
    
        if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
                print(f"Figure saved to {save_path}")
            
        return fig, ax

def create_box_plot(df,figPath,prs_types=['main','epi','cardio','combined'],phenotype_col='PHENOTYPE',orientation='vertical', dpi=300,figsize=(12, 6)):
    
    # Prepare dataframe
    df = df.copy()
    if 2 in df[phenotype_col].unique():
        df['Status'] = df[phenotype_col].map({2: 'Case', 1: 'Control'})
    else:
        df['Status'] = df[phenotype_col].map({1: 'Case', 0: 'Control'}) 
        
    
    # Create subplots
    n_plots = len(prs_types)
    n_cols = 2
    n_rows = (n_plots + 1) // 2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_plots > 1 else [axes]
    
    for idx, prs_type in enumerate(prs_types):
        ax = axes[idx]
        
        # Determine column
        if prs_type == 'combined':
            bin_col = 'combined_centile_bin'
        else:
            bin_col = f'{prs_type}_centile_bin'
            
        if bin_col not in df.columns:
            ax.text(0.5, 0.5, f'Data not available\nfor {prs_type}',
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=12)
            continue
        
        # Prepare data
        plot_data = df[[bin_col, 'Status']].dropna()
        plot_data.rename(columns={bin_col: 'Percentile_Bin'}, inplace=True)
        
        # Create boxplot
        color = COHORT_COLORS.get(prs_type, 'gray')
        label = labels.get(prs_type, f'PRS_{prs_type}')
        
        sns.boxplot(data=plot_data, x='Status', y='Percentile_Bin', ax=ax,
                   palette={'Control': '#95A5A6', 'Case': color},
                   showfliers=False, linewidth=1.5)
        
        # Add threshold line
        ax.axhline(y=800, color='black', linestyle=':', linewidth=1.5, 
                  alpha=0.5)
        
        # Statistics annotation
        control_median = plot_data[plot_data['Status'] == 'Control']['Percentile_Bin'].median()
        case_median = plot_data[plot_data['Status'] == 'Case']['Percentile_Bin'].median()
        
#       ax.text(0.03, 0.08, f'Control median: {control_median:.0f}\nCase median: {case_median:.0f}',
#              transform=ax.transAxes, fontsize=9, verticalalignment='top',
#              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        # Formatting
        ax.set_xlabel('Status', fontsize=11, fontweight='bold')
        ax.set_ylabel('PRS Percentile Bin', fontsize=11, fontweight='bold')
        ax.set_title(label, fontsize=13, fontweight='bold', pad=10)
        ax.grid(True, alpha=0.3, axis='y')
        
    # Hide unused subplots
    for idx in range(len(prs_types), len(axes)):
        axes[idx].set_visible(False)
        
    plt.suptitle('T2D PRS Distribution: Cases vs Controls by Model', 
                fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    

    plt.savefig(figPath, dpi=dpi, bbox_inches='tight')
    print(f"Subplot boxplot saved to {figPath}")
    
    return fig, axes

    
def plot_prevalence_with_ci(df, prs_types=['main', 'epi', 'cardio', 'combined'],
                                                        phenotype_col='PHENOTYPE',
                                                        figsize=(12, 6), style='errorbar', 
                                                        n_display_bins=10,
                                                        show_threshold=True, threshold_percentile=80,
                                                        save_path=None, dpi=300):
        """
        Create prevalence plot with confidence intervals for PRS calculations.
        
        Parameters:
        -----------
        df : DataFrame
                Data with centile_bin columns and phenotype
        prs_types : list
                PRS types to plot
        phenotype_col : str
                Column for case/control status
        figsize : tuple
                Figure size
        style : str
                'shaded' for filled CI regions or 'errorbar' for error bars
        n_display_bins : int
                Number of bins to display (default 100, groups 10 bines each from 1000)
        show_threshold : bool
                Whether to show high-risk threshold line
        threshold_percentile : int
                Percentile for high-risk threshold (default 80 for top 20%)
        save_path : str
                Path to save figure (optional)
        dpi : int
                Resolution for saved figure
        
        Returns:
        --------
        fig, ax : matplotlib figure and axis objects
        """
    
        # Set up plot style
        sns.set_style("whitegrid")

    
        # Calculate overall prevalence for reference line
        if 2 in df[phenotype_col].unique():
                overall_prevalence = (df[phenotype_col] - 1).mean()
        else:
                overall_prevalence = df[phenotype_col].mean()
            
        # Create figure
        fig, ax = plt.subplots(figsize=figsize,sharey=True)
    
        bin_stats_all = pd.DataFrame()
        # Plot each PRS type
        for prs_type in prs_types:
                try:
                        # Calculate bin prevalences with CIs
                        bin_stats = calculate_bin_prevalence_with_ci(
                                df, prs_type, phenotype_col=phenotype_col,n_display_bins=n_display_bins
                        )
                        bin_stats['prs_type'] = prs_type
                        bin_stats_all = pd.concat([bin_stats_all,bin_stats],ignore_index=True)
                    
                        color = COHORT_COLORS.get(prs_type, 'gray')
                        label = labels.get(prs_type, f'PRS_{prs_type}')
                    
                        if style == 'shaded':
                                # Plot line with shaded CI region
                                ax.plot(bin_stats['display_bin'], bin_stats['prevalence'], 
                                                color=color, linewidth=2, label=label, zorder=3)
                                ax.fill_between(bin_stats['display_bin'], 
                                                                bin_stats['ci_lower'], 
                                                                bin_stats['ci_upper'],
                                                                color=color, alpha=0.2, zorder=2)
                            
                        elif style == 'errorbar':
                                # Plot with error bars
                                yerr_lower = bin_stats['prevalence'] - bin_stats['ci_lower']
                                yerr_upper = bin_stats['ci_upper'] - bin_stats['prevalence']
                                ax.errorbar(bin_stats['display_bin'], bin_stats['prevalence'],
                                                        yerr=[yerr_lower, yerr_upper],
                                                        color=color, linewidth=2, label=label,
                                                        capsize=3, capthick=1, alpha=0.8, zorder=3)
                            
                except Exception as e:
                        print(f"Warning: Could not plot {prs_type}: {e}")
                        continue
            
        # Add overall prevalence reference line
        ax.axhline(y=overall_prevalence, color='gray', linestyle='--', 
                            linewidth=1.5, alpha=0.7, label='Overall prevalence', zorder=1)
    
        # Add high-risk threshold line
        if show_threshold:
                ax.axvline(x=threshold_percentile, color='black', linestyle=':', 
                                    linewidth=2, alpha=0.5, label=f'Top {100-threshold_percentile}% threshold',
                                    zorder=1)
            
        # Formatting
        ax.set_xlabel('PRS Percentile', fontsize=12, fontweight='bold')
        ax.set_ylabel('T2D Prevalence', fontsize=12, fontweight='bold')
        ax.set_title('T2D Prevalence by PRS Percentile with 95% Confidence Intervals', 
                                fontsize=14, fontweight='bold', pad=5)
    
        # Set x-axis limits and ticks
        ax.set_xlim(.8, 10.2)
        ax.set_xticks(range(1, 11))
    
        # Format y-axis as percentage
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1%}'))
    
        # Legend
        ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    
        # Grid
        ax.grid(True, alpha=0.3, zorder=0)
    
        plt.tight_layout()
    
        # Save if path provided
        if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
                print(f"Figure saved to {save_path}")
            
        return fig, ax, bin_stats_all


def plot_prevalence_comparison_subplots(df, prs_types=['main', 'epi', 'cardio', 'combined'],
                                                            phenotype_col='PHENOTYPE', n_display_bins=10,
                                                            style="shaded", dpi=300,sharey=True, ylim=None,
                                                            figsize=(14, 10), save_path=None):
        """
        Create separate subplot for each PRS type with prevalence and CIs.
        
        Parameters:
        -----------
        df : DataFrame
                Data with centile_bin columns
        prs_types : list
                PRS types to plot
        phenotype_col : str
                Column for case/control status
        bin_n_display_bins : int
                Bin grouping width into n bins
        figsize : tuple
                Figure size
        save_path : str
                Path to save figure
        
        Returns:
        --------
        fig, axes : matplotlib figure and axes objects
        """
    
        # Set up colors and labels

    
        # Calculate overall prevalence
        if 2 in df[phenotype_col].unique():
                overall_prevalence = (df[phenotype_col] - 1).mean()
        else:
                overall_prevalence = df[phenotype_col].mean()
            
        # Create subplots
        n_plots = len(prs_types)
        n_cols = 2
        n_rows = (n_plots + 1) // 2
    
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize,sharey=True)
        axes = axes.flatten() if n_plots > 1 else [axes]
    
        # Calculate y-axis limits from all data
        all_prevalences = []
        all_cis = []
        
        for prs_type in prs_types:
            try:
                bin_stats = calculate_bin_prevalence_with_ci(
                    df, prs_type, phenotype_col=phenotype_col, n_display_bins=n_display_bins
                )
                all_prevalences.extend(bin_stats['prevalence'].values)
                all_cis.extend(bin_stats['ci_lower'].values)
                all_cis.extend(bin_stats['ci_upper'].values)
            except:
                pass
                
        if all_cis:
            y_min = max(0, min(all_cis) - 0.02)
            y_max = min(1, max(all_cis) + 0.02)
        else:
            y_min, y_max = 0, 0.5
    
        for idx, prs_type in enumerate(prs_types):
                ax = axes[idx]
            
                try:
                        # Calculate prevalence with CIs
                        bin_stats = calculate_bin_prevalence_with_ci(
                                df, prs_type, phenotype_col=phenotype_col,
                                 n_display_bins=n_display_bins
                        )
                    
                        color = COHORT_COLORS.get(prs_type, 'gray')
                        label = labels.get(prs_type, f'PRS_{prs_type}')
                    
                        if style == 'shaded':
                                
                            # Plot line with shaded CI
                            ax.plot(bin_stats['display_bin'], bin_stats['prevalence'], 
                                            color=color, linewidth=2.5, label=label)
                            ax.fill_between(bin_stats['display_bin'], 
                                                            bin_stats['ci_lower'], 
                                                            bin_stats['ci_upper'],
                                                            color=color, alpha=0.3)
                        else:
                            # Plot with error bars
                            yerr_lower = bin_stats['prevalence'] - bin_stats['ci_lower']
                            yerr_upper = bin_stats['ci_upper'] - bin_stats['prevalence']
                            ax.errorbar(bin_stats['display_bin'], bin_stats['prevalence'],
                                                    yerr=[yerr_lower, yerr_upper],
                                                    color=color, linewidth=2, label=label,
                                                    capsize=3, capthick=1, alpha=0.8, zorder=3)
                    
                        # Reference line
                        ax.axhline(y=overall_prevalence, color='gray', linestyle='--', 
                                            linewidth=1, alpha=0.7, label='Overall prevalence')
                    
                        # Threshold line
                        ax.axvline(x=80, color='black', linestyle=':', 
                                            linewidth=1.5, alpha=0.5)
                    
                        # Formatting
                        # Only show y-label on leftmost plot
                        if idx == 0:
                            ax.set_ylabel('T2D Prevalence', fontsize=11, fontweight='bold')
                        ax.set_xlabel('PRS Percentile', fontsize=11, fontweight='bold')
                        ax.set_title(label, fontsize=13, fontweight='bold', pad=3)
                        ax.set_xlim(.8, 10.2)
                        ax.set_ylim(y_min, y_max)
                        ax.set_xticks(range(1, 11))  # Changed from range(0, 10, 1) to range(1, 11)
                        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1%}'))
                        if idx == 0:
                            ax.legend(loc='upper left', fontsize=9)
                        ax.grid(True, alpha=0.3)
                    
                except Exception as e:
                        print(f"Warning: Could not plot {prs_type}: {e}")
                        ax.text(0.5, 0.5, f'Data not available\nfor {prs_type}', 
                                        ha='center', va='center', transform=ax.transAxes)
                    
        # Hide unused subplots
        for idx in range(len(prs_types), len(axes)):
                axes[idx].set_visible(False)
            
        plt.suptitle('T2D Prevalence by PRS Type with 95% Confidence Intervals', 
                                fontsize=15, fontweight='bold', y=0.995)
        plt.tight_layout()
    
        if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight')
                print(f"Figure saved to {save_path}")
            
        return fig, axes


def calculate_bin_prevalence_with_ci(df, prs_type, phenotype_col='PHENOTYPE', n_display_bins=100,
                                     n_bins=1000, method='wilson'):
    """
    Calculate prevalence with confidence intervals for bins or bin groups.
    
    Parameters:
    -----------
    df : DataFrame
        Data with centile_bin and phenotype columns
    prs_type : str
        PRS type ('main', 'epi', 'cardio', 'combined')
    phenotype_col : str
        Column name for case/control status
    n_bins : int
        Total number of bins (default 1000)
    n_display_bins : int
        Number of bins to display default n_bins
    bin_width : int
        Width of bin groups for smoother plotting (default 10 bins)
    method : str
        CI method: 'wilson', 'binomial', or 'bootstrap'
    
    Returns:
    --------
    bin_stats : DataFrame
        Prevalence and CIs for each bin group
    """
    
    # Adjust phenotype if needed
    if 2 in df[phenotype_col].unique():
        phenotype = df[phenotype_col] - 1
    else:
        phenotype = df[phenotype_col]
        
    # Determine bin column name
    if prs_type == 'combined':
        bin_col = 'combined_centile_bin'
    else:
        bin_col = f'{prs_type}_centile_bin'
        
    if bin_col not in df.columns:
        raise ValueError(f"Column {bin_col} not found in dataframe")
        
    # Group bins for smoother plotting
#   bin_groups = range(1, n_bins + 1, bin_width)
    bin_width = n_bins // n_display_bins  # For 100 display bins: 1000/100 = 10 bins each
    
    
    results = []
    
    for display_bin in range(1, n_display_bins + 1):
        # Calculate range of original bins for this display bin
        bin_start = (display_bin - 1) * bin_width + 1
#       bin_end = min(bin_start + bin_width - 1, n_bins)
        bin_end = display_bin * bin_width
        
        # Get individuals in this bin range
        mask = (df[bin_col] >= bin_start) & (df[bin_col] <= bin_end)
        subset_phenotype = phenotype[mask]
        
        if len(subset_phenotype) == 0:
            continue
        
        # Calculate prevalence
        n_cases = subset_phenotype.sum()
        n_total = len(subset_phenotype)
        prevalence = n_cases / n_total
        
        # Calculate confidence interval
        if method == 'wilson':
            ci_lower, ci_upper = proportion_confint(n_cases, n_total, 
                                                   alpha=0.05, method='wilson')
        elif method == 'binomial':
            ci_lower, ci_upper = proportion_confint(n_cases, n_total, 
                                                   alpha=0.05, method='binom_test')
        elif method == 'bootstrap':
            # Bootstrap CI for this bin group
            if n_total > 10:  # Only bootstrap if enough samples
                bootstrap_prevs = []
                for _ in range(1000):
                    resampled = np.random.choice(subset_phenotype.values, 
                                               size=len(subset_phenotype), 
                                               replace=True)
                    bootstrap_prevs.append(resampled.mean())
                ci_lower = np.percentile(bootstrap_prevs, 2.5)
                ci_upper = np.percentile(bootstrap_prevs, 97.5)
            else:
                ci_lower, ci_upper = prevalence, prevalence
                
        # Calculate bin midpoint for plotting
#       bin_midpoint = (bin_start + bin_end) / 2
#       percentile = bin_midpoint / 10  # Convert to percentile (0-100)
        percentile = (display_bin - 0.5)  # Percentile 0.5 to 100.5
        
        
#       results.append({
#           'bin_start': bin_start,
#           'bin_end': bin_end,
#           'bin_midpoint': bin_midpoint,
#           'percentile': percentile,
#           'prevalence': prevalence,
#           'ci_lower': ci_lower,
#           'ci_upper': ci_upper,
#           'n_cases': int(n_cases),
#           'n_total': int(n_total)
#       })
#       
#       
#   return pd.DataFrame(results)
        
        results.append({
            'display_bin': display_bin,
            'bin_start': bin_start,
            'bin_end': bin_end,
            'percentile': percentile,
            'prevalence': prevalence,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_cases': int(n_cases),
            'n_total': int(n_total)
        })
        
    return pd.DataFrame(results)


def calculate_mean_prevalence_across_bins(df, phenotype_col='PHENOTYPE', prs_types=['main', 'epi', 'cardio'],
                                                                method='bootstrap', n_bootstrap=10000):
        """
        Calculate mean prevalence across bins for each PRS type with confidence intervals.
        
        Parameters:
        -----------
        df : DataFrame
                DataFrame with columns: IID, PHENOTYPE, {type}_centile_bin, {type}_high_risk, etc.
        phenotype_col : str
                Column name for case/control status
        prs_types : list
                List of PRS types to analyze (e.g., ['main', 'epi', 'cardio'])
        method : str
                'bootstrap' or 'normal' for CI calculation
        n_bootstrap : int
                Number of bootstrap iterations
        
        Returns:
        --------
        prevalence_stats : DataFrame
                Mean prevalence and CIs for each PRS type across all bins
        """
    
        # Adjust phenotype if coded as 2/1 instead of 1/0
        if 2 in df[phenotype_col].unique():
                phenotype = df[phenotype_col] - 1
        else:
                phenotype = df[phenotype_col]
            
        results = []
    
        for prs_type in prs_types:
                bin_col = f'{prs_type}_centile_bin'
            
                if bin_col not in df.columns:
                        print(f"Warning: {bin_col} not found, skipping {prs_type}")
                        continue
            
                # Calculate prevalence for each bin
                bin_prevalences = []
                for bin_num in range(1, 1001):  # Bins 1-1000
                        mask = df[bin_col] == bin_num
                        if mask.sum() > 0:
                                bin_prevalence = phenotype[mask].mean()
                                bin_prevalences.append(bin_prevalence)
                            
                # Convert to array
                case_rates = np.array(bin_prevalences)
            
                # Calculate mean prevalence across bins
                mean_prevalence = np.mean(case_rates)
            
                if method == 'bootstrap':
                        # Bootstrap confidence intervals
                        bootstrap_means = []
                        for _ in range(n_bootstrap):
                                resampled = np.random.choice(case_rates, size=len(case_rates), replace=True)
                                bootstrap_means.append(np.mean(resampled))
                            
                        ci_lower = np.percentile(bootstrap_means, 2.5)
                        ci_upper = np.percentile(bootstrap_means, 97.5)
                    
                elif method == 'normal':
                        # Normal approximation
                        std_error = stats.sem(case_rates)
                        ci_margin = std_error * stats.t.ppf(0.975, len(case_rates) - 1)
                        ci_lower = mean_prevalence - ci_margin
                        ci_upper = mean_prevalence + ci_margin
                    
                results.append({
                        'prs_type': f'PRS_{prs_type}',
                        'mean_prevalence': mean_prevalence,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper,
                        'n_bins': len(case_rates),
                        'method': method
                })
            
        # Add composite PRS
        composite_prevalences = []
        for bin_num in range(1, 1001):
                mask = df['combined_centile_bin'] == bin_num
                if mask.sum() > 0:
                        bin_prevalence = phenotype[mask].mean()
                        composite_prevalences.append(bin_prevalence)
                    
        case_rates_composite = np.array(composite_prevalences)
        mean_prevalence_composite = np.mean(case_rates_composite)
    
        if method == 'bootstrap':
                bootstrap_means = []
                for _ in range(n_bootstrap):
                        resampled = np.random.choice(case_rates_composite, size=len(case_rates_composite), replace=True)
                        bootstrap_means.append(np.mean(resampled))
                ci_lower = np.percentile(bootstrap_means, 2.5)
                ci_upper = np.percentile(bootstrap_means, 97.5)
        else:
                std_error = stats.sem(case_rates_composite)
                ci_margin = std_error * stats.t.ppf(0.975, len(case_rates_composite) - 1)
                ci_lower = mean_prevalence_composite - ci_margin
                ci_upper = mean_prevalence_composite + ci_margin
            
        results.append({
                'prs_type': 'PRS_comp',
                'mean_prevalence': mean_prevalence_composite,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'n_bins': len(case_rates_composite),
                'method': method
        })
    
        return pd.DataFrame(results)


def calculate_high_risk_prevalence(df, phenotype_col='PHENOTYPE',
                                prs_types=['main', 'epi', 'cardio']):
        """
        Calculate prevalence in high-risk groups (top 20%) with confidence intervals.
        
        Parameters:
        -----------
        df : DataFrame
                DataFrame with {type}_high_risk columns (0/1)
        phenotype_col : str
                Column name for case/control status
        prs_types : list
                List of PRS types to analyze
        
        Returns:
        --------
        high_risk_stats : DataFrame
                Prevalence statistics for high-risk groups
        """
        from statsmodels.stats.proportion import proportion_confint
    
        # Adjust phenotype if needed
        if 2 in df[phenotype_col].unique():
                phenotype = df[phenotype_col] - 1
        else:
                phenotype = df[phenotype_col]
            
        results = []
    
        for prs_type in prs_types:
                high_risk_col = f'{prs_type}_high_risk'
            
                if high_risk_col not in df.columns:
                        continue
            
                # Filter to high-risk individuals
                high_risk_mask = df[high_risk_col] == 1
                high_risk_phenotype = phenotype[high_risk_mask]
            
                if len(high_risk_phenotype) == 0:
                        continue
            
                # Calculate prevalence
                n_cases = high_risk_phenotype.sum()
                n_total = len(high_risk_phenotype)
                prevalence = n_cases / n_total
            
                # Wilson score interval
                ci_lower, ci_upper = proportion_confint(n_cases, n_total, alpha=0.05, method='wilson')
            
                results.append({
                        'prs_type': f'PRS_{prs_type}',
                        'prevalence': prevalence,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper,
                        'n_cases': int(n_cases),
                        'n_total': int(n_total),
                        'percentile': 'Top 20%'
                })
            
        # Add composite
        high_risk_mask = df['combined_high_risk'] == 1
        high_risk_phenotype = phenotype[high_risk_mask]
    
        n_cases = high_risk_phenotype.sum()
        n_total = len(high_risk_phenotype)
        prevalence = n_cases / n_total
        ci_lower, ci_upper = proportion_confint(n_cases, n_total, alpha=0.05, method='wilson')
    
        results.append({
                'prs_type': 'PRS_crmult',
                'prevalence': prevalence,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'n_cases': int(n_cases),
                'n_total': int(n_total),
                'percentile': 'Top 20%'
        })
    
        return pd.DataFrame(results)


def calculate_prevalence_by_percentile_ranges(df, phenotype_col='PHENOTYPE',
                                                            prs_types=['main', 'epi', 'cardio'],
                                                            ranges=None):
        """
        Calculate prevalence for custom percentile ranges.
        
        Parameters:
        -----------
        df : DataFrame
                DataFrame with centile_bin columns
        phenotype_col : str
                Column for case/control status
        prs_types : list
                PRS types to analyze
        ranges : list of dict
                Custom ranges, e.g., [{'name': 'Top 10%', 'min_bin': 900, 'max_bin': 1000}]
                If None, uses default ranges (top 5%, 10%, 20%)
        
        Returns:
        --------
        range_stats : DataFrame
                Prevalence for each range and PRS type
        """
        from statsmodels.stats.proportion import proportion_confint
    
        if ranges is None:
                ranges = [
                        {'name': 'Top 3%', 'min_bin': 970, 'max_bin': 1000},
                        {'name': 'Top 5%', 'min_bin': 950, 'max_bin': 1000},
                        {'name': 'Top 7%', 'min_bin': 930, 'max_bin': 1000},
                        {'name': 'Top 10%', 'min_bin': 900, 'max_bin': 1000},
                        {'name': 'Top 20%', 'min_bin': 800, 'max_bin': 1000}
                ]
            
        # Adjust phenotype
        if 2 in df[phenotype_col].unique():
                phenotype = df[phenotype_col] - 1
        else:
                phenotype = df[phenotype_col]
            
        results = []
    
        for prs_type in prs_types + ['combined']:
                if prs_type == 'combined':
                        bin_col = 'combined_centile_bin'
                        prs_label = 'PRS_combined'
                else:
                        bin_col = f'{prs_type}_centile_bin'
                        prs_label = f'PRS_{prs_type}'
                    
                if bin_col not in df.columns:
                        continue
            
                for range_dict in ranges:
                        mask = (df[bin_col] >= range_dict['min_bin']) & (df[bin_col] <= range_dict['max_bin'])
                        subset_phenotype = phenotype[mask]
                    
                        if len(subset_phenotype) == 0:
                                continue
                    
                        n_cases = subset_phenotype.sum()
                        n_total = len(subset_phenotype)
                        prevalence = n_cases / n_total
                    
                        ci_lower, ci_upper = proportion_confint(n_cases, n_total, 
                                                                                                        alpha=0.05, method='wilson')
                    
                        results.append({
                                'prs_type': prs_label,
                                'percentile_range': range_dict['name'],
                                'prevalence': prevalence,
                                'ci_lower': ci_lower,
                                'ci_upper': ci_upper,
                                'n_cases': int(n_cases),
                                'n_total': int(n_total)
                        })
                    
        return pd.DataFrame(results)


# Example usage:
if __name__ == "__main__":
    
    pheno_data = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi'
    holdoutPath = f'{pheno_data}/scores/combinedPRSGroups.holdout.csv'
    figPath = f'{pheno_data}/figures'
    
#   # Load your data
#   holdout_df = pd.read_csv(holdoutPath)
#
#   print("="*70)
#   print("MEAN PREVALENCE ACROSS ALL BINS (1-1000)")
#   print("="*70)
#
#   # Calculate mean prevalence across all bins
#   prevalence_stats = calculate_mean_prevalence_across_bins(
#           holdout_df, 
#           prs_types=['main', 'epi', 'cardio'],
#           method='bootstrap',
#           n_bootstrap=10000
#   )
#
#   print("\nMean prevalence across bins (Bootstrap 95% CI):")
#   print(prevalence_stats.to_string(index=False))
#
#   print("\n" + "="*70)
#   print("HIGH-RISK GROUP PREVALENCE (Top 20%, bins ≥800)")
#   print("="*70)
#
#   # Calculate prevalence in high-risk groups
#   high_risk_stats = calculate_high_risk_prevalence(
#           holdout_df,
#           prs_types=['main', 'epi', 'cardio']
#   )
#
#   print("\nPrevalence in high-risk groups (Wilson 95% CI):")
#   print(high_risk_stats.to_string(index=False))
#
#   print("\n" + "="*70)
#   print("PREVALENCE BY PERCENTILE RANGES")
#   print("="*70)
#
#   # Calculate prevalence for different percentile ranges
#   range_stats = calculate_prevalence_by_percentile_ranges(
#           holdout_df,
#           prs_types=['main', 'epi', 'cardio']
#   )
#
#   print("\nPrevalence by percentile range (Wilson 95% CI):")
#   print(range_stats.to_string(index=False))
#
#   # Create formatted output table
#   print("\n" + "="*70)
#   print("SUMMARY TABLE FOR MANUSCRIPT")
#   print("="*70)
#
#   summary = prevalence_stats.copy()
#   summary['formatted_ci'] = summary.apply(
#           lambda x: f"{x['mean_prevalence']:.3f} ({x['ci_lower']:.3f}-{x['ci_upper']:.3f})",
#           axis=1
#   )
#   print("\nMean prevalence across bins [95% CI]:")
#   print(summary[['prs_type', 'formatted_ci']].to_string(index=False))
#   
#
#       
#   print("Creating prevalence plots...")
#   
#
#   # Single plot with all PRS types
#   fig1, ax1, bin_stats = plot_prevalence_with_ci(
#           holdout_df,
#           prs_types=['main', 'epi', 'cardio', 'combined'],
#           phenotype_col='PHENOTYPE',
#           style='errorbar',
#           figsize=(12, 6),
#           dpi=300,
#           save_path=f'{figPath}/prevalence_all_prs_with_ci.png'
#   )
#   
#   bin_stats.to_csv(f'{pheno_data}/scores/combinedPRSGroups.prevalence.holdout.csv',index=False)
#   
#   # Separate subplots
#   fig2, axes2 = plot_prevalence_comparison_subplots(
#           holdout_df,
#           prs_types=['main', 'epi', 'cardio', 'combined'],
#           style="shaded",
#           phenotype_col='PHENOTYPE',
#           figsize=(14, 10),
#           save_path=f'{figPath}/prevalence_subplots_with_ci.png'
#   )
#
#   # Custom plot for just composite vs baseline
#   fig3, ax3, bin_stats = plot_prevalence_with_ci(
#           holdout_df,
#           prs_types=['main', 'combined'],
#           phenotype_col='PHENOTYPE',
#           style='errorbar',
#           figsize=(10, 6),
#           save_path=f'{figPath}/prevalence_composite_vs_baseline.png'
#   )
#   
#   fig4, ax4 = create_box_plot(holdout_df,f'{figPath}/box_plots_composite_vs_all.png',prs_types=['main','epi','cardio','combined'],phenotype_col='PHENOTYPE',figsize=(16, 12))
#   
#   # Separate subplots
#   fig5, ax5 = plot_prevalence_comparison_subplots(
#           holdout_df,
#           prs_types=['main','combined'],
#           style="shaded",
#           phenotype_col='PHENOTYPE',
#           figsize=(14, 10),
#           save_path=f'{figPath}/prevalence_subplots_with_ci_main_combined.png'
#   )
    
    or_df = pd.read_csv(f'{pheno_data}/scores/combinedORPRSGroups.holdout.csv')
    fig6,ax6 = plot_odds_ratios_comparison(or_df, percentiles=[1,3,5,10,20],prs_methods=['main','combined'], figsize=(12, 8),
                                        show_significance=True, reference_line=True,
                                        save_path=f'{figPath}/odds_ratio_main_combined.png', dpi=300)
    
    