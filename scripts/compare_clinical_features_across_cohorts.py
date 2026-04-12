#!/usr/bin/env python3
"""
compare_clinical_features_across_cohorts.py

Compares raw clinical/environmental features across PRS-derived risk cohorts
using holdout (preferred) or validation individuals.

Design decisions:
  - Uses RAW untransformed features for clinical interpretability
  - Handles missing data per-feature (listwise deletion per test, not per person)
  - Includes features with >5% missingness (excluded from modelling but valid here)
  - Non-parametric tests (Mann-Whitney U) given non-normal raw clinical data
  - FDR correction across features within each comparison
  - Compares: each high-risk cohort vs all other high-risk cohorts
              each low-risk cohort vs all other low-risk cohorts
              high-risk vs low-risk within each cohort
"""

import pandas as pd
import numpy as np
import os
import argparse
import warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")


# =============================================================================
#  STATISTICAL COMPARISON UTILITIES
# =============================================================================

def compare_two_groups(series_a, series_b, feature_name, label_a, label_b):
    """
    Compare two groups on a single continuous feature using Mann-Whitney U.
    Returns a dict with test results, effect size, and summary statistics.
    Uses only non-missing values for each group independently.
    """
    a = series_a.dropna()
    b = series_b.dropna()

    n_a = len(a)
    n_b = len(b)

    if n_a < 5 or n_b < 5:
        return None  # Too few observations for a meaningful test

    try:
        stat, pval = stats.mannwhitneyu(a, b, alternative='two-sided')
        # Rank-biserial correlation as effect size (r = 1 - 2U / n_a*n_b)
        effect_size = 1 - (2 * stat) / (n_a * n_b)
    except Exception:
        return None

    return {
        'feature': feature_name,
        'group_a': label_a,
        'group_b': label_b,
        'n_a': n_a,
        'n_b': n_b,
        'mean_a': a.mean(),
        'median_a': a.median(),
        'std_a': a.std(),
        'mean_b': b.mean(),
        'median_b': b.median(),
        'std_b': b.std(),
        'mean_diff': a.mean() - b.mean(),
        'median_diff': a.median() - b.median(),
        'mw_statistic': stat,
        'p_value': pval,
        'effect_size_r': effect_size,   # rank-biserial r: >0 = a higher, <0 = b higher
        'missing_a': series_a.isna().sum(),
        'missing_b': series_b.isna().sum(),
        'pct_missing_a': series_a.isna().mean() * 100,
        'pct_missing_b': series_b.isna().mean() * 100,
    }


def run_comparisons(data, feature_cols, group_col, groups, comparison_type):
    """
    Run all pairwise comparisons between cohorts across all features.

    comparison_type: 'each_vs_others' (one cohort vs all remaining combined)
                     'pairwise'        (all unique pairs)
    """
    rows = []

    for feature in feature_cols:
        if comparison_type == 'each_vs_others':
            for g in groups:
                mask_g = data[group_col] == g
                mask_other = (data[group_col] != g) & data[group_col].isin(groups)
                result = compare_two_groups(
                    data.loc[mask_g, feature],
                    data.loc[mask_other, feature],
                    feature, g, f'other_{comparison_type}'
                )
                if result:
                    result['cohort'] = g
                    result['comparison_type'] = comparison_type
                    rows.append(result)

        elif comparison_type == 'pairwise':
            for i, g1 in enumerate(groups):
                for g2 in groups[i + 1:]:
                    mask_g1 = data[group_col] == g1
                    mask_g2 = data[group_col] == g2
                    result = compare_two_groups(
                        data.loc[mask_g1, feature],
                        data.loc[mask_g2, feature],
                        feature, g1, g2
                    )
                    if result:
                        result['cohort'] = f'{g1}_vs_{g2}'
                        result['comparison_type'] = comparison_type
                        rows.append(result)

    return pd.DataFrame(rows)


def apply_fdr(df, alpha=0.05):
    """Apply Benjamini-Hochberg FDR correction per cohort comparison."""
    if df.empty:
        return df

    corrected = []
    for cohort, grp in df.groupby('cohort'):
        pvals = grp['p_value'].values
        _, pvals_corrected, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
        grp = grp.copy()
        grp['p_value_fdr'] = pvals_corrected
        grp['significant_fdr'] = pvals_corrected < alpha
        corrected.append(grp)

    return pd.concat(corrected, ignore_index=True)


# =============================================================================
#  MISSINGNESS REPORT
# =============================================================================

def summarise_missingness(data, feature_cols):
    """Report missingness per feature across the dataset."""
    rows = []
    for f in feature_cols:
        n_total = len(data)
        n_missing = data[f].isna().sum()
        rows.append({
            'feature': f,
            'n_total': n_total,
            'n_missing': n_missing,
            'n_present': n_total - n_missing,
            'pct_missing': n_missing / n_total * 100,
            'excluded_from_modelling': (n_missing / n_total) > 0.05
        })
    return pd.DataFrame(rows).sort_values('pct_missing', ascending=False)


# =============================================================================
#  VISUALISATION
# =============================================================================

def plot_significant_features(data, results_df, feature_cols, group_col,
                               groups, fig_path, tag, top_n=20):
    """
    Boxplots for the top N most significant features (by FDR p-value)
    in the each_vs_others comparison.
    """
    if results_df.empty or 'p_value_fdr' not in results_df.columns:
        return

    sig = results_df[results_df['significant_fdr']].copy()
    if sig.empty:
        print(f"  No significant features at FDR<0.05 for {tag}")
        return

    top_features = (sig.groupby('feature')['p_value_fdr']
                    .min()
                    .sort_values()
                    .head(top_n)
                    .index.tolist())

    n_features = len(top_features)
    ncols = 4
    nrows = int(np.ceil(n_features / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = axes.flatten() if n_features > 1 else [axes]

    plot_data = data[data[group_col].isin(groups)][[group_col] + top_features].copy()

    for ax, feature in zip(axes, top_features):
        # Collect per-cohort data, drop NaN per feature
        cohort_data = [plot_data.loc[plot_data[group_col] == g, feature].dropna().values
                       for g in groups]
        ax.boxplot(cohort_data, labels=groups, patch_artist=True,
                   medianprops=dict(color='black', linewidth=2))
        ax.set_title(feature, fontsize=9, fontweight='bold')
        ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=7)
        ax.set_ylabel('Raw value', fontsize=7)

        # Mark significance
        min_p = sig[sig['feature'] == feature]['p_value_fdr'].min()
        ax.set_xlabel(f'FDR p={min_p:.3e}', fontsize=7)

    # Hide unused axes
    for ax in axes[n_features:]:
        ax.set_visible(False)

    plt.suptitle(f'Top {n_features} Significant Features — {tag}', fontsize=11, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, f'significant_features_{tag}.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved boxplot → significant_features_{tag}.png")


def plot_missingness_heatmap(missingness_df, fig_path):
    """Bar chart of missingness rates across features."""
    fig, ax = plt.subplots(figsize=(12, max(4, len(missingness_df) * 0.25)))
    colours = ['#d62728' if x else '#1f77b4' for x in missingness_df['excluded_from_modelling']]
    ax.barh(missingness_df['feature'], missingness_df['pct_missing'], color=colours)
    ax.axvline(5, color='red', linestyle='--', linewidth=1, label='5% threshold')
    ax.set_xlabel('% Missing')
    ax.set_title('Feature Missingness (red = excluded from modelling)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'feature_missingness.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved missingness plot → feature_missingness.png")


def plot_effect_size_heatmap(results_df, fig_path, tag):
    """Heatmap of effect sizes (rank-biserial r) for significant features."""
    if results_df.empty or 'significant_fdr' not in results_df.columns:
        return

    sig_features = results_df[results_df['significant_fdr']]['feature'].unique()
    if len(sig_features) == 0:
        return

    pivot = results_df[results_df['feature'].isin(sig_features)].pivot_table(
        index='feature', columns='cohort', values='effect_size_r', aggfunc='mean')

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 1.5),
                                    max(4, len(pivot) * 0.4)))
    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.2f', annot_kws={'size': 7},
                linewidths=0.5, cbar_kws={'label': 'Rank-biserial r'})
    ax.set_title(f'Effect Size Heatmap — {tag}\n(positive = focal cohort higher)')
    ax.set_xlabel('Cohort')
    ax.set_ylabel('Feature')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, f'effect_size_heatmap_{tag}.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved effect size heatmap → effect_size_heatmap_{tag}.png")


# =============================================================================
#  MAIN
# =============================================================================

def main(pheno_data, raw_features_file, holdout_ids_file=None,
         validation_ids_file=None, use_set='holdout',
         alpha=0.05, min_effect_size=0.1,use_all=False):

    scoresPath = f'{pheno_data}/scores'
    figPath = f'{pheno_data}/figures/cohortClinicalComparison'
    outputPath = f'{scoresPath}/cohortClinicalComparison'
    os.makedirs(figPath, exist_ok=True)
    os.makedirs(outputPath, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Load PRS bin assignments (cohort labels for all individuals)
    # -------------------------------------------------------------------------
    if use_set == 'holdout':
      combined_prs = pd.read_csv(f'{scoresPath}/combinedPRSGroups.holdout.filtered.csv')
    else:
      combined_prs = pd.read_csv(f'{scoresPath}/combinedPRSGroups.filtered.csv')
    bin_cols = [c for c in combined_prs.columns if c.startswith('bin_') and 'combined' not in c]
    if not use_all:
      bin_cols = [c for c in bin_cols if 'all' not in c]
      
    cohort_names = [c.replace('bin_', '') for c in bin_cols]
    print(f"Cohorts found: {cohort_names}")

    combined_prs = combined_prs[['IID', 'PHENOTYPE'] + bin_cols].set_index('IID')

    # -------------------------------------------------------------------------
    # 2. Determine which individuals to use
    # -------------------------------------------------------------------------
    if use_set == 'holdout' and holdout_ids_file:
        id_df = pd.read_csv(holdout_ids_file)
        id_col = [c for c in id_df.columns if 'IID' in c or 'id' in c.lower()][0]
        use_ids = id_df[id_col].tolist()
        print(f"Using holdout set: {len(use_ids)} individuals")
    elif use_set == 'validation' and validation_ids_file:
        id_df = pd.read_csv(validation_ids_file)
        id_col = [c for c in id_df.columns if 'IID' in c or 'id' in c.lower()][0]
        use_ids = id_df[id_col].tolist()
        print(f"Using validation set: {len(use_ids)} individuals")
    else:
        # Fall back: use all individuals
        use_ids = combined_prs.index.tolist()
        print(f"No ID file provided — using all {len(use_ids)} individuals")

    prs_subset = combined_prs.loc[combined_prs.index.isin(use_ids)].copy()
    print(f"Individuals with PRS data in selected set: {len(prs_subset)}")

    # -------------------------------------------------------------------------
    # 3. Load raw clinical features
    # -------------------------------------------------------------------------
    raw_df = pd.read_csv(raw_features_file)
    try:
        raw_df.rename(columns={'Participant ID': 'IID'}, inplace=True)
    except Exception:
        pass
    raw_df.set_index('IID', inplace=True)

    # Restrict to individuals in selected set
    raw_df = raw_df.loc[raw_df.index.isin(use_ids)]
    print(f"Individuals with raw clinical data in selected set: {len(raw_df)}")

    # Identify feature columns — drop any ID/metadata columns
    exclude_cols = ['PHENOTYPE', 'SEX', 'age'] + bin_cols + cohort_names
    feature_cols = [c for c in raw_df.columns if c not in exclude_cols]
    print(f"Clinical features to compare: {len(feature_cols)}")

    # -------------------------------------------------------------------------
    # 4. Merge PRS bins with raw features
    # -------------------------------------------------------------------------
    merged = prs_subset.join(raw_df[feature_cols], how='inner')
    print(f"Individuals after merging PRS + clinical data: {len(merged)}")

    # -------------------------------------------------------------------------
    # 5. Assign cohort membership for high-risk, low-risk, and exclusive high-risk
    #    High/LowRisk: a person can appear in multiple cohorts (overlapping)
    #    ExclusiveHighRisk: a person must be in the top bin for exactly ONE cohort
    # -------------------------------------------------------------------------
    # Detect bin scale (deciles 1-10 vs quantiles 1-1000) and derive thresholds
    bin_max = merged[bin_cols].max().max()
    high_thresh = 0.8 * bin_max   # top 20 %
    low_thresh  = 0.2 * bin_max   # bottom 20 %

    for col in bin_cols:
        cohort = col.replace('bin_', '')
        merged[f'highrisk_{cohort}'] = (merged[col] > high_thresh).astype(int)
        merged[f'lowrisk_{cohort}']  = (merged[col] < low_thresh).astype(int)

    # Exclusive high-risk: above threshold in THIS cohort only
    for col in bin_cols:
        cohort     = col.replace('bin_', '')
        other_cols = [c for c in bin_cols if c != col]
        is_high_this  = merged[col] > high_thresh
        is_high_other = (
            merged[other_cols].gt(high_thresh).any(axis=1)
            if other_cols else pd.Series(False, index=merged.index)
        )
        merged[f'exclusive_{cohort}'] = (
            is_high_this & ~is_high_other
        ).astype(int)

    # -------------------------------------------------------------------------
    # 6. Missingness report
    # -------------------------------------------------------------------------
    print("\nCalculating missingness...")
    missingness = summarise_missingness(merged, feature_cols)
    missingness.to_csv(os.path.join(outputPath, 'feature_missingness.csv'), index=False)
    plot_missingness_heatmap(missingness, figPath)

    n_high_missing = (missingness['pct_missing'] > 5).sum()
    print(f"  {n_high_missing}/{len(feature_cols)} features have >5% missing "
          f"(included here but flagged)")

    # -------------------------------------------------------------------------
    # 7. Build comparison datasets: one row per individual per cohort membership
    #    Strategy: for each risk tier (high/low), create a long-format dataframe
    #    with a 'cohort' column indicating which cohort that individual belongs to
    # -------------------------------------------------------------------------
    all_results = []

    for risk_tier, prefix, phenotype_filter in [
        ('HighRisk',          'highrisk',  2),
        ('LowRisk',           'lowrisk',   1),
        ('ExclusiveHighRisk', 'exclusive', 2),
    ]:
        print(f"\n{'='*60}")
        print(f"  Comparing {risk_tier} cohorts")
        print(f"{'='*60}")

        # For each cohort, get the individuals in that risk tier
        # Restrict to the correct phenotype (cases for high, controls for low)
        tier_frames = []
        cohort_ns = {}

        for cohort in cohort_names:
            mask = (merged[f'{prefix}_{cohort}'] == 1)
            if phenotype_filter is not None:
                mask &= (merged['PHENOTYPE'] == phenotype_filter)
            cohort_data = merged[mask][feature_cols].copy()
            cohort_data['cohort'] = cohort
            cohort_ns[cohort] = len(cohort_data)
            tier_frames.append(cohort_data)

        print(f"  Cohort sizes: { {k: v for k, v in cohort_ns.items()} }")

        tier_df = pd.concat(tier_frames, ignore_index=False)

        # Each-vs-others comparison
        results = run_comparisons(
            tier_df, feature_cols, group_col='cohort',
            groups=cohort_names, comparison_type='each_vs_others'
        )

        if results.empty:
            print(f"  No results generated for {risk_tier}")
            continue

        results = apply_fdr(results, alpha=alpha)

        # Flag meaningful effect sizes
        results['meaningful_effect'] = (
            results['significant_fdr'] &
            (results['effect_size_r'].abs() >= min_effect_size)
        )

        results['risk_tier'] = risk_tier

        # Merge in missingness info
        results = results.merge(
            missingness[['feature', 'pct_missing', 'excluded_from_modelling']],
            on='feature', how='left'
        )

        all_results.append(results)

        # Summary of findings
        n_sig = results['significant_fdr'].sum()
        n_meaningful = results['meaningful_effect'].sum()
        print(f"  Significant features (FDR<{alpha}): {n_sig}")
        print(f"  Significant + meaningful effect (|r|>={min_effect_size}): {n_meaningful}")

        # Highlight epi_product specifically
        epi_prod_results = results[
            (results['cohort'] == 'epi_product') &
            results['significant_fdr']
        ].sort_values('p_value_fdr')

        if len(epi_prod_results) > 0:
            print(f"\n  epi_product significant features ({len(epi_prod_results)}):")
            print(epi_prod_results[['feature', 'mean_a', 'mean_b', 'effect_size_r',
                                     'p_value_fdr']].to_string(index=False))
        else:
            print(f"\n  epi_product: no significant clinical features found")

        # Plots
        tag = f"{risk_tier}_{use_set}"
        plot_significant_features(tier_df, results, feature_cols,
                                   'cohort', cohort_names, figPath, tag)
        plot_effect_size_heatmap(results, figPath, tag)

        # Save per-tier results
        results.sort_values(['cohort', 'p_value_fdr']).to_csv(
            os.path.join(outputPath, f'clinical_comparison_{risk_tier}_{use_set}.csv'),
            index=False
        )
        print(f"  Saved results → clinical_comparison_{risk_tier}_{use_set}.csv")

    # -------------------------------------------------------------------------
    # 8. Combined output, priority logic, and summary table
    # -------------------------------------------------------------------------
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(
            os.path.join(outputPath, f'clinical_comparison_all_{use_set}.csv'),
            index=False
        )

        # ── Priority logic: ExclusiveHighRisk supersedes HighRisk ──────────
        # When a clinical feature is significant (FDR < alpha, |r| >= min_effect)
        # in BOTH HighRisk and ExclusiveHighRisk for the same cohort, keep only
        # the ExclusiveHighRisk row (higher specificity).  If ExclusiveHighRisk
        # is not significant for that feature × cohort, retain the HighRisk row.
        hr  = combined[combined['risk_tier'] == 'HighRisk'].copy()
        exc = combined[combined['risk_tier'] == 'ExclusiveHighRisk'].copy()
        lr  = combined[combined['risk_tier'] == 'LowRisk'].copy()

        exc_sig = exc[exc['meaningful_effect']].set_index(['feature', 'cohort'])

        def _superseded(row):
            return (row['feature'], row['cohort']) in exc_sig.index

        hr['exclusive_supersedes'] = hr.apply(_superseded, axis=1)
        exc['exclusive_supersedes'] = False
        lr['exclusive_supersedes']  = False

        # Priority output: ExclusiveHighRisk rows where significant; HighRisk
        # rows only where ExclusiveHighRisk is not significant for that pair.
        hr_kept = hr[~hr['exclusive_supersedes']]
        priority = pd.concat([hr_kept, exc, lr], ignore_index=True)
        priority.sort_values(['cohort', 'risk_tier', 'p_value_fdr'], inplace=True)
        priority.to_csv(
            os.path.join(outputPath, f'clinical_comparison_priority_{use_set}.csv'),
            index=False
        )
        print(f"\n  Priority output ({len(priority)} rows): HighRisk rows superseded "
              f"by ExclusiveHighRisk = {hr['exclusive_supersedes'].sum()}")
        print(f"  Saved → clinical_comparison_priority_{use_set}.csv")

        # Concise summary: features significant for epi_product in either tier
        epi_summary = priority[
            (priority['cohort'] == 'epi_product') &
            priority['meaningful_effect']
        ][['risk_tier', 'feature', 'n_a', 'mean_a', 'median_a',
           'mean_b', 'median_b', 'effect_size_r',
           'p_value_fdr', 'pct_missing']].sort_values(['risk_tier', 'p_value_fdr'])

        epi_summary.to_csv(
            os.path.join(outputPath, f'epi_product_clinical_summary_{use_set}.csv'),
            index=False
        )

        print(f"\n{'='*60}")
        print(f"SUMMARY")
        print(f"{'='*60}")
        print(f"Total tests run       : {len(combined)}")
        print(f"Significant (FDR<{alpha}) : {combined['significant_fdr'].sum()}")
        print(f"Meaningful effect     : {combined['meaningful_effect'].sum()}")
        print(f"\nepi_product meaningful clinical features (priority): {len(epi_summary)}")
        if not epi_summary.empty:
            print(epi_summary[['risk_tier', 'feature', 'effect_size_r',
                                'p_value_fdr']].to_string(index=False))

        print(f"\nAll outputs saved to: {outputPath}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compare raw clinical features across PRS risk cohorts")

    parser.add_argument("--pheno_data", help="Path to pheno results directory")
    parser.add_argument("--raw_features_file", help="Path to raw (untransformed) clinical features CSV")
    parser.add_argument("--holdout_ids_file", default=None,
                        help="CSV with IIDs of holdout individuals (preferred)")
    parser.add_argument("--validation_ids_file", default=None,
                        help="CSV with IIDs of validation individuals")
    parser.add_argument("--use_set", default='holdout', choices=['holdout', 'validation', 'all'],
                        help="Which individual set to use (default: holdout)")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="FDR significance threshold (default: 0.05)")
    parser.add_argument("--min_effect_size", type=float, default=0.1,
                        help="Minimum rank-biserial r to flag as meaningful (default: 0.1)")
    parser.add_argument("--use_all", type=bool, default=False,
                        help="use_all flag if using validation set")

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    raw_features_file = args.raw_features_file or os.environ.get("RAW_FEATURES_FILE")

    if not pheno_data:
        raise ValueError("Provide --pheno_data or set PHENO_DATA env var")
    if not raw_features_file:
        raise ValueError("Provide --raw_features_file or set RAW_FEATURES_FILE env var")

    print(f"Pheno data path  : {pheno_data}")
    print(f"Raw features file: {raw_features_file}")
    print(f"Using set        : {args.use_set}")
    print(f"FDR alpha        : {args.alpha}")
    print(f"Min effect size  : {args.min_effect_size}")

    main(
        pheno_data=pheno_data,
        raw_features_file=raw_features_file,
        holdout_ids_file=args.holdout_ids_file,
        validation_ids_file=args.validation_ids_file,
        use_set=args.use_set,
        alpha=args.alpha,
        min_effect_size=args.min_effect_size,
        use_all=args.use_all
    )
  