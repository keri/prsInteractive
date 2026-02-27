#!/usr/bin/env python3
"""
run_cohort_analysis_pipeline.py

Orchestrates the full per-cohort analysis pipeline:
  Step 1 — calculate_top_features_in_cohort:
            Trains multi-output RF on high/low risk cohorts, extracts SHAP-based
            genomic/interaction features predictive of each cohort.

  Step 2 — compare_clinical_features_across_cohorts:
            Compares raw clinical/environmental features across PRS risk cohorts
            using Mann-Whitney U and FDR correction.

  Step 3 — compile_results:
            Joins outputs from both steps per cohort into unified summary tables
            and a cross-cohort comparison report. Cohorts with signal only in one
            source (e.g. epi_product: clinical signal but no genomic features)
            are explicitly flagged.

Outputs (all under {pheno_data}/scores/cohortSummary/):
  cohort_summary_{cohort}.csv          — full feature list for each cohort
  pipeline_overview.csv                — one row per cohort, high-level counts
  cohorts_with_signal_discordance.csv  — cohorts where genomic/clinical signals disagree
  figures/                             — combined plots per cohort
"""

import pandas as pd
import numpy as np
import os
import argparse
import warnings
import traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# Import both analysis modules
import calculate_top_features_in_cohort as genomic_module
import compare_clinical_features_across_cohorts as clinical_module


# =============================================================================
#  STEP 3: COMPILE RESULTS
# =============================================================================

def load_genomic_results(scores_path, filter_strategy):
    """Load SHAP-based feature importance results from step 1."""
    path = os.path.join(
        scores_path,
        'importantCohortScores',
        f'importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv'
    )
    if not os.path.exists(path):
        print(f"  [WARN] Genomic results file not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    df['source'] = 'genomic_shap'
    # Normalise column name: feature_set -> cohort for consistent merging
    if 'feature_set' in df.columns and 'cohort' not in df.columns:
        df.rename(columns={'feature_set': 'cohort'}, inplace=True)
    return df


def load_clinical_results(scores_path, use_set):
    """Load clinical comparison results from step 2."""
    path = os.path.join(
        scores_path,
        'cohortClinicalComparison',
        f'clinical_comparison_all_{use_set}.csv'
    )
    if not os.path.exists(path):
        print(f"  [WARN] Clinical results file not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    df['source'] = 'clinical_comparison'
    return df


def load_auc_results(scores_path):
    """Load model AUC performance table from step 1."""
    path = os.path.join(
        scores_path,
        'importantCohortScores',
        'performanceMetricsAcrossCohortsAndTrainingData.csv'
    )
    if not os.path.exists(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def build_cohort_summary(cohort, genomic_df, clinical_df, auc_df):
    """
    Build a unified summary for a single cohort combining genomic and clinical
    signal. Returns a dict of DataFrames:
      genomic   — SHAP features for this cohort
      clinical  — significant clinical features for this cohort
      overview  — single-row summary dict
    """
    # --- Genomic features for this cohort ---
    if not genomic_df.empty:
        cohort_genomic = genomic_df[genomic_df['cohort'] == cohort].copy()
    else:
        cohort_genomic = pd.DataFrame()

    # --- Clinical features for this cohort (significant only) ---
    if not clinical_df.empty:
        cohort_clinical = clinical_df[
            (clinical_df['cohort'] == cohort) &
            (clinical_df['significant_fdr'] == True)
        ].copy()
        cohort_clinical_all = clinical_df[clinical_df['cohort'] == cohort].copy()
    else:
        cohort_clinical = pd.DataFrame()
        cohort_clinical_all = pd.DataFrame()

    # --- AUC for this cohort ---
    auc_vals = {}
    if not auc_df.empty and 'index' in auc_df.columns:
        cohort_auc_rows = auc_df[auc_df['index'] == cohort]
        if not cohort_auc_rows.empty and 'auc' in cohort_auc_rows.columns:
            for _, row in cohort_auc_rows.iterrows():
                key = f"auc_{row.get('features_used','?')}_{row.get('training_data_used','?')}"
                auc_vals[key] = row['auc']

    # --- Overview row ---
    n_genomic_tier1 = len(cohort_genomic[cohort_genomic.get('specificity_tier', pd.Series(dtype=int)) == 1]) \
        if 'specificity_tier' in cohort_genomic.columns else 0
    n_genomic_total = len(cohort_genomic)
    n_clinical_highrisk = len(cohort_clinical[cohort_clinical['risk_tier'] == 'HighRisk']) \
        if 'risk_tier' in cohort_clinical.columns else 0
    n_clinical_lowrisk = len(cohort_clinical[cohort_clinical['risk_tier'] == 'LowRisk']) \
        if 'risk_tier' in cohort_clinical.columns else 0

    has_genomic_signal = n_genomic_total > 0
    has_clinical_signal = (n_clinical_highrisk + n_clinical_lowrisk) > 0
    signal_type = _classify_signal(has_genomic_signal, has_clinical_signal)

    overview = {
        'cohort': cohort,
        'n_genomic_features_total': n_genomic_total,
        'n_genomic_features_tier1': n_genomic_tier1,
        'n_clinical_features_highrisk': n_clinical_highrisk,
        'n_clinical_features_lowrisk': n_clinical_lowrisk,
        'has_genomic_signal': has_genomic_signal,
        'has_clinical_signal': has_clinical_signal,
        'signal_type': signal_type,
        **auc_vals
    }

    return {
        'genomic': cohort_genomic,
        'clinical_significant': cohort_clinical,
        'clinical_all': cohort_clinical_all,
        'overview': overview
    }


def _classify_signal(has_genomic, has_clinical):
    """Label the signal pattern for a cohort."""
    if has_genomic and has_clinical:
        return 'both'
    elif has_genomic and not has_clinical:
        return 'genomic_only'
    elif not has_genomic and has_clinical:
        return 'clinical_only'
    else:
        return 'no_signal'


# =============================================================================
#  PLOTS
# =============================================================================

def plot_cohort_overview(overview_df, fig_path):
    """Bar chart comparing genomic and clinical feature counts per cohort."""
    if overview_df.empty:
        return

    cohorts = overview_df['cohort'].tolist()
    x = np.arange(len(cohorts))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, len(cohorts) * 1.5), 5))

    bars1 = ax.bar(x - width, overview_df['n_genomic_features_total'],
                   width, label='Genomic features (all tiers)', color='#1f77b4', alpha=0.8)
    bars2 = ax.bar(x, overview_df['n_genomic_features_tier1'],
                   width, label='Genomic features (tier 1 only)', color='#aec7e8', alpha=0.8)
    bars3 = ax.bar(x + width, overview_df['n_clinical_features_highrisk'],
                   width, label='Clinical features (high risk)', color='#d62728', alpha=0.8)

    ax.set_xlabel('Cohort')
    ax.set_ylabel('Number of significant features')
    ax.set_title('Genomic vs Clinical Feature Signal per Cohort')
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, rotation=30, ha='right')
    ax.legend()

    # Annotate signal type
    signal_colours = {'both': 'green', 'genomic_only': 'blue',
                      'clinical_only': 'orange', 'no_signal': 'red'}
    for i, (_, row) in enumerate(overview_df.iterrows()):
        colour = signal_colours.get(row['signal_type'], 'grey')
        ax.annotate(row['signal_type'].replace('_', '\n'),
                    xy=(x[i], 0), xytext=(x[i], -0.8),
                    ha='center', fontsize=7, color=colour)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'cohort_feature_signal_overview.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved → cohort_feature_signal_overview.png")


def plot_auc_comparison(overview_df, fig_path):
    """Plot AUC values per cohort if available."""
    auc_cols = [c for c in overview_df.columns if c.startswith('auc_')]
    if not auc_cols:
        return

    melted = overview_df[['cohort'] + auc_cols].melt(
        id_vars='cohort', var_name='model', value_name='auc')
    melted.dropna(subset=['auc'], inplace=True)
    if melted.empty:
        return

    fig, ax = plt.subplots(figsize=(max(8, len(overview_df) * 1.5), 5))
    sns.barplot(data=melted, x='cohort', y='auc', hue='model', ax=ax, alpha=0.8)
    ax.axhline(0.5, color='red', linestyle='--', linewidth=1, label='Random (0.5)')
    ax.set_ylim(0.4, 1.0)
    ax.set_title('Model AUC per Cohort and Feature Set')
    ax.set_xlabel('Cohort')
    ax.set_ylabel('AUC')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'cohort_auc_comparison.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved → cohort_auc_comparison.png")


def plot_signal_discordance(discordance_df, fig_path):
    """Heatmap showing which cohorts have genomic vs clinical signal."""
    if discordance_df.empty:
        return

    signal_map = {'both': 2, 'genomic_only': 1, 'clinical_only': -1, 'no_signal': 0}
    plot_df = discordance_df[['cohort', 'signal_type']].copy()
    plot_df['signal_code'] = plot_df['signal_type'].map(signal_map)
    plot_df = plot_df.set_index('cohort')

    fig, ax = plt.subplots(figsize=(4, max(3, len(plot_df) * 0.5)))
    cmap = sns.diverging_palette(20, 145, n=5, as_cmap=True)
    sns.heatmap(plot_df[['signal_code']], annot=plot_df[['signal_type']].values,
                fmt='', cmap=cmap, center=0, ax=ax,
                cbar_kws={'label': 'Signal type'},
                linewidths=0.5)
    ax.set_title('Signal Source per Cohort')
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'signal_discordance_heatmap.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved → signal_discordance_heatmap.png")


def plot_top_genomic_features_per_cohort(genomic_df, fig_path, top_n=15):
    """Horizontal bar chart of top genomic features per cohort coloured by tier."""
    if genomic_df.empty:
        return

    tier_colours = {1: '#d62728', 2: '#ff7f0e', 3: '#1f77b4'}
    cohorts = genomic_df['cohort'].unique()

    for cohort in cohorts:
        cdf = genomic_df[genomic_df['cohort'] == cohort].copy()
        if cdf.empty:
            continue

        sort_col = 'matching_z_score' if 'matching_z_score' in cdf.columns else \
                   'odds_ratio' if 'odds_ratio' in cdf.columns else None
        if sort_col:
            cdf = cdf.sort_values(sort_col, ascending=False).head(top_n)

        colours = [tier_colours.get(t, 'grey')
                   for t in cdf.get('specificity_tier', pd.Series([3] * len(cdf)))]

        fig, ax = plt.subplots(figsize=(8, max(3, len(cdf) * 0.4)))
        bars = ax.barh(cdf['feature'], cdf[sort_col] if sort_col else range(len(cdf)),
                       color=colours, alpha=0.85)
        ax.set_xlabel(sort_col or 'rank')
        ax.set_title(f'Top Genomic Features — {cohort}')
        ax.invert_yaxis()

        # Tier legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=c, label=f'Tier {t}')
                           for t, c in tier_colours.items()]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

        plt.tight_layout()
        plt.savefig(os.path.join(fig_path, f'genomic_features_{cohort}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()

    print(f"  Saved genomic feature plots for {len(cohorts)} cohorts")


def plot_clinical_effect_sizes_per_cohort(clinical_df, fig_path):
    """Dot plot of effect sizes for significant clinical features per cohort."""
    if clinical_df.empty or 'significant_fdr' not in clinical_df.columns:
        return

    sig = clinical_df[clinical_df['significant_fdr'] == True].copy()
    if sig.empty:
        return

    cohorts = sig['cohort'].unique()

    for cohort in cohorts:
        cdf = sig[sig['cohort'] == cohort].copy()
        if cdf.empty:
            continue

        # Separate high and low risk
        for risk_tier in ['HighRisk', 'LowRisk']:
            tdf = cdf[cdf['risk_tier'] == risk_tier].copy() if 'risk_tier' in cdf.columns else cdf
            if tdf.empty:
                continue

            tdf = tdf.sort_values('effect_size_r', ascending=True)

            fig, ax = plt.subplots(figsize=(7, max(3, len(tdf) * 0.35)))
            colours = ['#d62728' if r > 0 else '#1f77b4' for r in tdf['effect_size_r']]
            ax.barh(tdf['feature'], tdf['effect_size_r'], color=colours, alpha=0.8)
            ax.axvline(0, color='black', linewidth=0.8)
            ax.set_xlabel('Rank-biserial r (positive = focal cohort higher)')
            ax.set_title(f'Clinical Feature Effect Sizes — {cohort} [{risk_tier}]')

            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, f'clinical_effects_{cohort}_{risk_tier}.png'),
                        dpi=150, bbox_inches='tight')
            plt.close()

    print(f"  Saved clinical effect size plots for {len(cohorts)} cohorts")


# =============================================================================
#  COMPILE MAIN
# =============================================================================

def compile_results(pheno_data, filter_strategy, use_set, output_path, fig_path):
    """
    Load outputs from both steps, build per-cohort summaries,
    and write compiled reports.
    """
    scores_path = f'{pheno_data}/scores'

    print("\n" + "="*60)
    print("STEP 3: COMPILING RESULTS")
    print("="*60)

    genomic_df = load_genomic_results(scores_path, filter_strategy)
    clinical_df = load_clinical_results(scores_path, use_set)
    auc_df = load_auc_results(scores_path)

    print(f"  Genomic features loaded : {len(genomic_df)} rows")
    print(f"  Clinical features loaded: {len(clinical_df)} rows")

    if genomic_df.empty and clinical_df.empty:
        print("  [ERROR] No results found from either step. Cannot compile.")
        return

    # Get full cohort list from whichever source has data
    cohorts = set()
    if not genomic_df.empty and 'cohort' in genomic_df.columns:
        cohorts |= set(genomic_df['cohort'].unique())
    if not clinical_df.empty and 'cohort' in clinical_df.columns:
        cohorts |= set(clinical_df['cohort'].unique())
    cohorts = sorted(cohorts)
    print(f"  Cohorts to summarise: {cohorts}")

    overview_rows = []
    all_genomic = []
    all_clinical = []

    for cohort in cohorts:
        print(f"\n  Summarising cohort: {cohort}")
        summary = build_cohort_summary(cohort, genomic_df, clinical_df, auc_df)

        overview_rows.append(summary['overview'])

        # Annotate with cohort for stacking
        if not summary['genomic'].empty:
            all_genomic.append(summary['genomic'])

        if not summary['clinical_significant'].empty:
            all_clinical.append(summary['clinical_significant'])

        # Write per-cohort file combining both sources
        cohort_genomic = summary['genomic'].copy()
        cohort_clinical = summary['clinical_significant'].copy()

        cohort_genomic['analysis'] = 'genomic_shap'
        cohort_clinical['analysis'] = 'clinical_comparison'

        # Align columns for stacking — use common columns where possible
        shared_cols = ['cohort', 'feature', 'analysis', 'training_data_used'] \
            if 'training_data_used' in cohort_genomic.columns else ['cohort', 'feature', 'analysis']

        genomic_out = cohort_genomic.reindex(columns=shared_cols + [
            c for c in cohort_genomic.columns if c not in shared_cols])
        clinical_out = cohort_clinical.reindex(columns=shared_cols + [
            c for c in cohort_clinical.columns if c not in shared_cols])

        cohort_combined = pd.concat([genomic_out, clinical_out], ignore_index=True)
        cohort_file = os.path.join(output_path, f'cohort_summary_{cohort}.csv')
        cohort_combined.to_csv(cohort_file, index=False)

        sig = summary['overview']['signal_type']
        n_g = summary['overview']['n_genomic_features_total']
        n_c = summary['overview']['n_clinical_features_highrisk']
        print(f"    signal_type={sig}  genomic={n_g}  clinical_highrisk={n_c}")

    # ----- Overview table -----
    overview_df = pd.DataFrame(overview_rows)
    overview_df.to_csv(os.path.join(output_path, 'pipeline_overview.csv'), index=False)
    print(f"\n  Saved pipeline_overview.csv")

    # ----- Discordance table -----
    discordance_df = overview_df[overview_df['signal_type'] != 'both'].copy()
    discordance_df.to_csv(
        os.path.join(output_path, 'cohorts_with_signal_discordance.csv'), index=False)
    print(f"  Saved cohorts_with_signal_discordance.csv ({len(discordance_df)} cohorts)")

    # ----- All genomic features stacked -----
    if all_genomic:
        stacked_genomic = pd.concat(all_genomic, ignore_index=True)
        stacked_genomic.to_csv(
            os.path.join(output_path, 'all_cohorts_genomic_features.csv'), index=False)
        print(f"  Saved all_cohorts_genomic_features.csv ({len(stacked_genomic)} rows)")

    # ----- All clinical features stacked -----
    if all_clinical:
        stacked_clinical = pd.concat(all_clinical, ignore_index=True)
        stacked_clinical.to_csv(
            os.path.join(output_path, 'all_cohorts_clinical_features.csv'), index=False)
        print(f"  Saved all_cohorts_clinical_features.csv ({len(stacked_clinical)} rows)")

    # ----- Feature overlap: which features appear in both analyses -----
    if all_genomic and all_clinical:
        genomic_features = set(pd.concat(all_genomic)['feature'].dropna().unique())
        clinical_features = set(pd.concat(all_clinical)['feature'].dropna().unique())
        overlap = genomic_features & clinical_features

        if overlap:
            print(f"\n  Features appearing in BOTH genomic and clinical analyses: {len(overlap)}")
            overlap_df = pd.DataFrame({'feature': sorted(overlap)})
            overlap_df['in_genomic'] = True
            overlap_df['in_clinical'] = True
            overlap_df.to_csv(
                os.path.join(output_path, 'genomic_clinical_feature_overlap.csv'), index=False)
        else:
            print("\n  No feature name overlap between genomic and clinical analyses "
                  "(expected — different feature spaces)")

    # ----- Plots -----
    print("\n  Generating summary plots...")
    plot_cohort_overview(overview_df, fig_path)
    plot_auc_comparison(overview_df, fig_path)
    plot_signal_discordance(discordance_df, fig_path)

    if all_genomic:
        plot_top_genomic_features_per_cohort(
            pd.concat(all_genomic, ignore_index=True), fig_path)
    if all_clinical:
        plot_clinical_effect_sizes_per_cohort(
            pd.concat(all_clinical, ignore_index=True), fig_path)

    # ----- Print final summary -----
    print("\n" + "="*60)
    print("PIPELINE COMPLETE — SUMMARY")
    print("="*60)
    print(overview_df[['cohort', 'signal_type', 'n_genomic_features_total',
                        'n_clinical_features_highrisk',
                        'n_clinical_features_lowrisk']].to_string(index=False))
    print(f"\nAll outputs saved to: {output_path}")


# =============================================================================
#  ORCHESTRATOR
# =============================================================================

def run_pipeline(pheno_data, feature_scores_file, raw_features_file,
                 holdout_ids_file=None, validation_ids_file=None,
                 use_set='holdout', threshold=1.99,
                 filter_strategy='tiered', majority_threshold=0.6,
                 delta_threshold=0.5, alpha=0.05, min_effect_size=0.1,
                 skip_genomic=False, skip_clinical=False,use_all=False):

    output_path = os.path.join(pheno_data, 'scores', 'cohortSummary')
    fig_path = os.path.join(pheno_data, 'figures', 'cohortSummary')
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(fig_path, exist_ok=True)

    step_status = {}

    # -------------------------------------------------------------------------
    # STEP 1: Genomic SHAP feature importance
    # -------------------------------------------------------------------------
    if not skip_genomic:
        print("\n" + "="*60)
        print("STEP 1: GENOMIC FEATURE IMPORTANCE (calculate_top_features_in_cohort)")
        print("="*60)
        try:
            genomic_module.main(
                phenoData=pheno_data,
                feature_scores_file=feature_scores_file,
                threshold=threshold,
                filter_strategy=filter_strategy,
                majority_threshold=majority_threshold,
                delta_threshold=delta_threshold
            )
            step_status['genomic'] = 'success'
            print("  Step 1 completed successfully")
        except Exception as e:
            step_status['genomic'] = f'failed: {e}'
            print(f"  [ERROR] Step 1 failed: {e}")
            traceback.print_exc()
    else:
        step_status['genomic'] = 'skipped'
        print("\nStep 1 skipped (--skip_genomic)")

    # -------------------------------------------------------------------------
    # STEP 2: Clinical feature comparison
    # -------------------------------------------------------------------------
    if not skip_clinical:
        print("\n" + "="*60)
        print("STEP 2: CLINICAL FEATURE COMPARISON (compare_clinical_features_across_cohorts)")
        print("="*60)
        try:
            clinical_module.main(
                pheno_data=pheno_data,
                raw_features_file=raw_features_file,
                holdout_ids_file=holdout_ids_file,
                validation_ids_file=validation_ids_file,
                use_set=use_set,
                alpha=alpha,
                min_effect_size=min_effect_size,
                use_all=use_all
            )
            step_status['clinical'] = 'success'
            print("  Step 2 completed successfully")
        except Exception as e:
            step_status['clinical'] = f'failed: {e}'
            print(f"  [ERROR] Step 2 failed: {e}")
            traceback.print_exc()
    else:
        step_status['clinical'] = 'skipped'
        print("\nStep 2 skipped (--skip_clinical)")

    # -------------------------------------------------------------------------
    # STEP 3: Compile results (runs regardless of step 1/2 success
    #         as long as output files exist from a prior run)
    # -------------------------------------------------------------------------
    compile_results(
        pheno_data=pheno_data,
        filter_strategy=filter_strategy,
        use_set=use_set,
        output_path=output_path,
        fig_path=fig_path
    )

    # Final status report
    print("\n" + "="*60)
    print("STEP STATUS")
    print("="*60)
    for step, status in step_status.items():
        print(f"  {step:12s}: {status}")
    print(f"  compile     : completed")


# =============================================================================
#  CLI
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Full cohort analysis pipeline: genomic features + clinical comparison + compilation",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required
    parser.add_argument("--pheno_data",
                        help="Path to pheno results directory")
    parser.add_argument("--feature_scores_file_filtered",
                        help="Path to feature scores file with beta coefficients (for step 1)")
    parser.add_argument("--raw_features_file",
                        help="Path to raw (untransformed) clinical features CSV (for step 2)")

    # Individual set selection
    parser.add_argument("--holdout_ids_file", default=None,
                        help="CSV with IIDs of holdout individuals")
    parser.add_argument("--validation_ids_file", default=None,
                        help="CSV with IIDs of validation individuals")
    parser.add_argument("--use_set", default='holdout',
                        choices=['holdout', 'validation', 'all'],
                        help="Which individual set to use for clinical comparison (default: holdout)")

    # Genomic step parameters
    parser.add_argument("--threshold", type=float, default=1.99,
                        help="SHAP z-score threshold per cohort (default: 1.99)")
    parser.add_argument("--filter_strategy", default='tiered',
                        choices=['strict', 'majority', 'differential', 'tiered'],
                        help="Cross-cohort filtering strategy (default: tiered)")
    parser.add_argument("--majority_threshold", type=float, default=0.6,
                        help="Majority strategy: fraction of other cohorts that must be negative (default: 0.6)")
    parser.add_argument("--delta_threshold", type=float, default=0.5,
                        help="Differential strategy: min delta z above mean of others (default: 0.5)")

    # Clinical step parameters
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="FDR alpha for clinical comparisons (default: 0.05)")
    parser.add_argument("--min_effect_size", type=float, default=0.1,
                        help="Min rank-biserial r to flag as meaningful (default: 0.1)")
    parser.add_argument("--use_all", type=bool, default=False,
                        help="include model with all features in analysis, default to False")

    # Pipeline control
    parser.add_argument("--skip_genomic", action='store_true',
                        help="Skip step 1 and use existing genomic results (e.g. recompile only)")
    parser.add_argument("--skip_clinical", action='store_true',
                        help="Skip step 2 and use existing clinical results (e.g. recompile only)")
      

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    feature_scores_file = args.feature_scores_file_filtered or os.environ.get("FEATURE_SCORES_FILE_FILTERED")
    raw_features_file = args.raw_features_file or os.environ.get("RAW_FEATURES_FILE")

    if not pheno_data:
        raise ValueError("Provide --pheno_data or PHENO_DATA env var")
    if not feature_scores_file and not args.skip_genomic:
        raise ValueError("Provide --feature_scores_file_filtered or FEATURE_SCORES_FILE_FILTERED env var")
    if not raw_features_file and not args.skip_clinical:
        raise ValueError("Provide --raw_features_file or RAW_FEATURES_FILE env var")

    print("="*60)
    print("COHORT ANALYSIS PIPELINE")
    print("="*60)
    print(f"  pheno_data          : {pheno_data}")
    print(f"  feature_scores_file : {feature_scores_file}")
    print(f"  raw_features_file   : {raw_features_file}")
    print(f"  use_set             : {args.use_set}")
    print(f"  use_all             : {args.use_all}")
    print(f"  filter_strategy     : {args.filter_strategy}")
    print(f"  threshold           : {args.threshold}")
    print(f"  majority_threshold  : {args.majority_threshold}")
    print(f"  delta_threshold     : {args.delta_threshold}")
    print(f"  alpha               : {args.alpha}")
    print(f"  min_effect_size     : {args.min_effect_size}")
    print(f"  skip_genomic        : {args.skip_genomic}")
    print(f"  skip_clinical       : {args.skip_clinical}")

    run_pipeline(
        pheno_data=pheno_data,
        feature_scores_file=feature_scores_file,
        raw_features_file=raw_features_file,
        holdout_ids_file=args.holdout_ids_file,
        validation_ids_file=args.validation_ids_file,
        use_set=args.use_set,
        threshold=args.threshold,
        filter_strategy=args.filter_strategy,
        majority_threshold=args.majority_threshold,
        delta_threshold=args.delta_threshold,
        alpha=args.alpha,
        min_effect_size=args.min_effect_size,
        skip_genomic=args.skip_genomic,
        skip_clinical=args.skip_clinical,
        use_all=args.use_all
    )
  