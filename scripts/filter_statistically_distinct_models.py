"""
filter_statically_distinct_models.py

Identifies statistically distinct PRS models, filters data to those models,
runs the full PRS plot suite, and generates the metrics comparison figure.

Usage
-----
    from filter_statically_distinct_models import filter_statically_distinct_models

    models_to_keep, filtered_data = filter_statically_distinct_models('/data/pheno_xyz')
"""

import pandas as pd
import argparse
import os
from helper.prs_model_filter import filter_redundant_models, quick_filter_models, apply_model_filter
from helper.draw_plots import create_all_prs_plots
from plot_prs_metrics_comparison import plot_prs_metrics_comparison, plot_performance_complementarity


def filter_statistically_distinct_models(pheno_data, use_all=None):
    """
    Identify statistically distinct PRS models, filter data, and produce all figures.

    Parameters
    ----------
    pheno_data : str
        Root phenotype data directory.  Expected layout:
            {pheno_data}/scores/
                McNemarStatsTestsAcrossPRSCalculations_refactored.csv
                model_recall_precision_improvement.csv
                combinedPRSGroups.csv
            {pheno_data}/figures/modelComparisons/   ← all outputs land here

    use_all : bool or None
        If truthy, include 'all v combined' in the pairwise filtering.
        If falsy / None (default), that comparison is excluded.

    Returns
    -------
    tuple[list, pd.DataFrame]
        (models_to_keep, filtered_data)
    """
    scores_path = f'{pheno_data}/scores'
    fig_path    = f'{pheno_data}/figures/modelComparisons'

    mcnemar_path = f'{scores_path}/McNemarStatsTestsAcrossPRSCalculations_refactored.csv'
    perf_path    = f'{scores_path}/model_recall_precision_improvement.csv'

    # ── 1. Decide which comparisons to include ─────────────────────────────
    if use_all:
        exclude_comparisons = None
    else:
        exclude_comparisons = ['all v combined']

    # ── 2. Statistical filtering — identify distinct models ────────────────
    results = quick_filter_models(
        mcnemar_path,
        exclude_comparisons=exclude_comparisons,
        output_path=f'{scores_path}/pairwise_filtering_decisions.csv'
    )

    models_to_keep = results['models_to_keep']
    print(f"Keep these models: {models_to_keep}")

    # ── 3. Filter PRS data to kept models only ─────────────────────────────
    prs_data      = pd.read_csv(f'{scores_path}/combinedPRSGroups.csv')
    filtered_data = apply_model_filter(prs_data, models_to_keep)
    filtered_data.to_csv(f'{scores_path}/combinedPRSGroups.filtered.csv', index=False)

    # ── 4. PRS correlation / high-risk plots ──────────────────────────────
    # models_to_keep already excludes prs_all if use_all=False (via filtering logic)
    # but we pass use_all to plots for consistency and potential future use
    create_all_prs_plots(
        combined_prs       = filtered_data,
        models_to_keep     = models_to_keep,
        output_path        = fig_path,
        threshold          = 8,
        priority_order     = ['cardio', 'epi', 'epi+main', 'main'],
        create_pairwise    = False,
        create_comprehensive = True,
        create_matrix      = True,
        include_all        = bool(use_all),  # match the filtering flag
        comprehensive_axes = None,
    )

    # ── 5. Metrics comparison figure ───────────────────────────────────────
    # Passes models_to_keep so Panels B and C show only the distinct models.
    # The heatmap still includes one representative row per redundant model
    # so the exclusion decision remains visible.
    plot_prs_metrics_comparison(
        mcnemar_path        = mcnemar_path,
        perf_path           = perf_path,
        output_path         = fig_path,
        models_to_keep      = models_to_keep,
        exclude_comparisons = exclude_comparisons,
        verbose             = True,
    )

    # ── 6. Performance & complementarity figure ────────────────────────────
    # Two-panel figure: Panel A = % cases identified (bars) + AUC twin axis;
    # Panel B = precision / recall grouped bars.
    # For a cross-phenotype comparison (e.g. T2D + Celiac side-by-side), call
    # plot_performance_complementarity() directly with multi_phenotype=; see
    # its docstring in plot_prs_metrics_comparison.py for full details.
    plot_performance_complementarity(
        mcnemar_path    = mcnemar_path,
        perf_path       = perf_path,
        output_path     = fig_path,
        models_to_keep  = models_to_keep,
        phenotype_label = os.path.basename(os.path.dirname(pheno_data)),
        verbose         = True,
    )

    # ── 6. (Optional) QQ plot — uncomment when ready ──────────────────────
    # from create_qq_plot_groups import create_qq_plot_groups
    # create_qq_plot_groups(filtered_data, fig_path)

    return models_to_keep, filtered_data


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="calculating filtering statistically distinct models and creating figures....")
    parser.add_argument("--pheno_data", help="Path to the input pheno folder")
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
#   pheno='type2Diabetes'
#   pheno_data = f'/Users/kerimulterer/prsInteractive/results/{pheno}/combinedAnalysis'
    
    models_to_keep, filtered_data = filter_statistically_distinct_models(pheno_data, use_all=None)
    