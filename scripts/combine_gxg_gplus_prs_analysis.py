#!/usr/bin/env python3
"""
combine_gxg_gplus_prs_analysis.py

Merges PRS results from GxG (productEpi) and G+G (summedEpi) analyses into a
single unified dataset, then replicates the complete downstream pipeline that
runs on each analysis separately:

    Step 1 – Merge & stage combined data
        • Loads combinedPRSGroups.csv from both analysis directories.
        • ``main`` is shared → kept from productEpi only (one copy).
        • Non-main models are renamed with ``_product`` / ``_summed`` suffixes.
        • Stages a combined prsScores/ directory with correctly-named files so
          calculate_top_features_in_cohort can glob them as expected.
          Source file patterns:
            test    – ``{model}.{n_features}.mixed.prs.csv``
            holdout – ``{model}.{n_features}.holdout.mixed.prs.csv``
        • Merges McNemar stats and performance CSVs from both sources.

    Step 2 – Scale & bin combined PRS columns
        • Replicates the scaling / binning logic in combine_prs.py.
        • ``prs_all`` models (``all_product`` / ``all_summed``) are excluded by
          default.  Pass ``use_all=True`` (or ``--use_all`` on the CLI)
          to include them.  Changing this flag requires re-running Steps 2–6.
        • ``combined_prs`` carry-overs from source analyses are stripped in
          Step 1 and recalculated fresh in Step 4 (assign_holdout_risk_bin).
        • Writes combinedPRSGroups.csv, combinedPrevalencePRSGroups.csv, and
          combinedORPRSGroups.csv to the combined scores directory.

    Step 3 – filter_statistically_distinct_models
        • McNemar tests and precision/recall calculations run on ALL models
          present in ``combinedPRSGroups.csv`` (including ``prs_all`` models).
        • Generates ``combinedPRSGroups.filtered.csv`` with statistically
          distinct models only.
        • ``use_all`` flag controls whether "all v combined" comparisons are
          included in the filtering logic. Default: ``False`` (exclude "all v
          combined" from filtering, so ``prs_all`` models won't be kept even
          if statistically distinct). Set ``True`` to allow ``prs_all`` models
          in the final filtered set.
        • NOTE: "combined" refers to the single-analysis meta-combination scores
          from productEpi/summedEpi holdout sets. These are filtered out during
          CSV merge in Step 2 since they don't exist in combinedAnalysis
          validation set.

    Step 4 – assign_holdout_risk_bin
        • Reads combinedPRSGroups.filtered.csv (test set, post-filter) to build
          1,000-bin risk statistics per scaled PRS column.
        • Assigns each holdout individual to a risk bin using training-derived
          bin edges, computes combined_centile_bin / combined_high_risk, and
          writes combinedPRSGroups.holdout.filtered.csv.

    Step 5 – calculate_top_features_in_cohort
        • SHAP-based feature importance across the distinct combined models.

    Step 6 – clinical_measure_performance
        • AUC / NRI analysis of the combined PRS against clinical measures.

Directory layout produced
-------------------------
    {pheno_results}/combinedAnalysis/
        scores/
            combinedPRSGroups.csv
            combinedPRSGroups.filtered.csv
            combinedPrevalencePRSGroups.csv
            combinedORPRSGroups.csv
            pairwise_filtering_decisions.csv
            McNemarStatsTestsAcrossPRSCalculations_refactored.csv   (merged)
            model_recall_precision_improvement.csv                  (merged)
            prsScores/          ← symlinks or copies, correctly renamed
            importantCohortScores/
        figures/
            modelComparisons/
            importantCohortFeatures/
            clinicalPerformance/

Usage
-----
    python combine_gxg_gplus_prs_analysis.py \\
        --pheno_results /path/to/{pheno}/results \\
        --pheno         type2Diabetes \\
        --results_path  /path/to/results \\
        --feature_scores_file /path/to/featureScoresReducedFinalModel.filtered.csv \\
        [--use_all] \\
        [--threshold 1.99] \\
        [--steps 1,2,3,4,5,6]     # run all steps by default

Notes
-----
* The script imports from the same helper modules used by the individual
  scripts (helper.draw_plots, helper.data_wrangling, etc.).  Make sure the
  Python path includes the project root so these imports resolve.
* ``filter_statistically_distinct_models``, ``assign_holdout_risk_bin``,
  ``calculate_top_features_in_cohort``, and ``clinical_measure_performance``
  are called via their ``main()`` functions —
  no subprocess overhead.
"""

import argparse
import glob
import os
import re
import shutil
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", message="Bad circle positioning")

# ── Local project imports ────────────────────────────────────────────────────
# Adjust sys.path if the project root is not already on it.
# e.g. sys.path.insert(0, '/path/to/project/root')

from combine_prs import (
    bin_prs,
    calculate_prevalance,
    scale_holdout_data_manually,
)
from filter_statistically_distinct_models import filter_statistically_distinct_models
from assign_holdout_risk_bin import main as run_holdout_bins
from run_cohort_analysis_pipeline import run_pipeline
from clinical_measure_performance import run_clinical_analysis as run_clinical_performance
from calculate_prs_stats import calculate_prs_stats

from helper.data_wrangling import (
    scale_data,
    calculate_training_statistics,
    calculate_odds_ratio_for_prs,
)
from helper.draw_plots import create_optimized_prevalence_plot

# ── Naming conventions ───────────────────────────────────────────────────────
PRODUCT_EPI  = "productEpi"
SUMMED_EPI   = "summedEpi"
COMBINED_DIR = "combinedAnalysis"

# ``main`` comes from both analyses but is identical → keep one copy only.
SHARED_MODELS = {"main"}

# Non-main models that exist in each analysis and will be suffixed.
METHOD_MODELS = {"epi", "epi+main", "cardio", "all"}

# The "all" model is optional — excluded by default (use_all=False).
# It represents a combined feature set and may not be present in every run.
ALL_MODEL = "all"

PRODUCT_SUFFIX = "_product"
SUMMED_SUFFIX  = "_summed"

# Columns produced by assign_holdout_risk_bin that must NOT be carried over
# from a prior single-analysis run into the combined merge.  These columns are
# created ONLY for the holdout set in Step 4 and should never appear in
# combinedPRSGroups.csv (the validation/test set), but we guard against carry-over
# defensively in case a user manually edited the CSV or ran a custom script.
#
# NOTE: prs_combined / scaled_prs_combined are NOT included here — they exist
# only in the holdout set and are used for scatter plots (e.g.,
# 'all_high_risk_cases_main_vs_epi.png'), never in the validation CSV.
# Columns created by assign_holdout_risk_bin that must not carry over
COMBINED_CARRYOVER_COLS = re.compile(
    r"^(bin_combined_centile_bin|combined_high_risk|bin_combined"
    r"any_high_risk|n_case_bins|n_valid_prs|prop_case_bins)$"
)

# ── Visual identity for the combined model set ───────────────────────────────
COHORT_COLORS_COMBINED: dict = {
    "main":              "#E69F00",   # orange  (shared)
    # productEpi models
    "epi_product":       "#56B4E9",   # sky blue
    "epi+main_product":  "#CC79A7",   # pinkish-purple
    "cardio_product":    "#009E73",   # bluish-green
    "all_product":       "#F0E442",   # yellow
    # summedEpi models
    "epi_summed":        "#0072B2",   # deep blue
    "epi+main_summed":   "#D55E00",   # vermillion
    "cardio_summed":     "#117733",   # dark green
    "all_summed":        "#AA4499",   # purple
}

COHORT_MARKERS_COMBINED: dict = {
    "main":              "o",
    "epi_product":       "s",
    "epi+main_product":  "p",
    "cardio_product":    "^",
    "all_product":       "x",
    "epi_summed":        "D",
    "epi+main_summed":   "h",
    "cardio_summed":     "v",
    "all_summed":        "+",
}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_dirs(*paths: str) -> None:
    """Create directories (including parents) silently."""
    for p in paths:
        os.makedirs(p, exist_ok=True)


def _model_rename_map(suffix: str) -> dict:
    """
    Return a column-rename dict for non-shared PRS columns.

    Applied to DataFrame columns (not filenames); the ``{n_features}``
    segment lives only in the filename, not in column names.

    e.g. suffix='_product' maps
        'prs_epi'        -> 'prs_epi_product'
        'scaled_prs_epi' -> 'scaled_prs_epi_product'
        'bin_epi'        -> 'bin_epi_product'
    """
    rename = {}
    for m in METHOD_MODELS:
        for prefix in ("prs_", "scaled_prs_", "bin_"):
            old = f"{prefix}{m}"
            new = f"{prefix}{m}{suffix}"
            rename[old] = new
    return rename


def _load_combined_prs(scores_path: str) -> pd.DataFrame:
    """Load combinedPRSGroups.csv from a scores directory."""
    path = os.path.join(scores_path, "combinedPRSGroups.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"combinedPRSGroups.csv not found at: {path}")
    return pd.read_csv(path)


def _suffix_model_name(model: str, suffix: str) -> str:
    """Apply suffix to a bare model name if it is not a shared model."""
    return model if model in SHARED_MODELS else f"{model}{suffix}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 – Merge & stage combined data
# ─────────────────────────────────────────────────────────────────────────────

def merge_prs_data(
    product_scores: str,
    summed_scores: str,
    combined_scores: str,
    product_prs_models: list[str] = None,
    summed_prs_models: list[str] = None,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Merge PRS group data from productEpi and summedEpi.

    ``main`` is taken from productEpi only; all other model columns are
    renamed with ``_product`` / ``_summed`` suffixes before merging.

    Returns
    -------
    merged_df : pd.DataFrame   — full combined dataset (unscaled)
    prs_list  : list[str]      — bare model names present in merged_df
    """
    print("\n[Step 1] Merging PRS data from productEpi and summedEpi ...")
    
    df_product = _load_combined_prs(product_scores)
      
    df_summed = _load_combined_prs(summed_scores)
      

    # ── Identify the PRS columns in each source ───────────────────────────
    def _get_raw_prs_cols(df: pd.DataFrame, prs_list: list[str]) -> list[str]:
      if prs_list is None:
        # Return all prs columns if no filter list provided
        return [c for c in df.columns 
            if c.startswith("prs_") and "scaled" not in c]
      return [c for c in df.columns 
          if c.startswith("prs_") and "scaled" not in c 
          and c.replace('prs_','') in prs_list]

    raw_product = _get_raw_prs_cols(df_product,product_prs_models)
    raw_summed  = _get_raw_prs_cols(df_summed,summed_prs_models)
    
    df_product = df_product[['IID','PHENOTYPE'] + raw_product]
    df_summed = df_summed[['IID','PHENOTYPE'] + raw_summed]
  
    # ── Drop columns we will recompute or must not carry over ───────────────
    # • scaled_prs_* and bin_* are recomputed in Step 2
    # • COMBINED_CARRYOVER_COLS (combined_centile_bin, prs_combined, etc.) were
    #   produced by a prior single-analysis assign_holdout_risk_bin or combine_prs
    #   run and must be recalculated for the combined dataset in Steps 2 and 4.
    drop_pattern = re.compile(r"^(scaled_prs_|bin_)")
    def _drop_recomputed(df: pd.DataFrame) -> pd.DataFrame:
        keep = [
            c for c in df.columns
            if not drop_pattern.match(c) and not COMBINED_CARRYOVER_COLS.match(c)
        ]
        dropped = set(df.columns) - set(keep)
        if dropped:
            print(f"    Dropped carry-over columns: {sorted(dropped)}")
        return df[keep]

    df_product = _drop_recomputed(df_product)
    df_summed  = _drop_recomputed(df_summed)

    # ── Rename non-main columns in each source ────────────────────────────
    df_product.rename(columns=_model_rename_map(PRODUCT_SUFFIX), inplace=True)
    df_summed.rename(columns=_model_rename_map(SUMMED_SUFFIX),  inplace=True)

    # ── Identify which summed columns are truly non-shared (drop main dup) ─
    summed_non_shared_cols = [
        c for c in df_summed.columns
        if c.startswith("prs_") and not any(c == f"prs_{m}" for m in SHARED_MODELS)
    ]
    merge_cols = ["IID", "PHENOTYPE"] + summed_non_shared_cols

    # ── Merge on IID + PHENOTYPE ─────────────────────────────────────────
    merged = df_product.merge(
        df_summed[merge_cols],
        on=["IID", "PHENOTYPE"],
        how="inner",
        suffixes=("", "_dup"),   # guard against any unexpected duplicates
    )
    # Drop any accidental duplicate columns
    merged = merged[[c for c in merged.columns if not c.endswith("_dup")]]

    # ── Build the final bare model list in display order ─────────────────
    prs_list = (
        [m for m in SHARED_MODELS]
        + [f"{m}{PRODUCT_SUFFIX}" for m in sorted(METHOD_MODELS)]
        + [f"{m}{SUMMED_SUFFIX}"  for m in sorted(METHOD_MODELS)]
    )
    # Keep only models that actually ended up with a column
    prs_list = [m for m in prs_list if f"prs_{m}" in merged.columns]

    n_shared = len(merged)
    print(f"  productEpi rows : {len(df_product):,}")
    print(f"  summedEpi rows  : {len(df_summed):,}")
    print(f"  inner-join rows : {n_shared:,}")
    print(f"  models retained : {prs_list}")

    return merged, prs_list

def _read_models_to_keep_from_filtering(scores_path: str) -> list[str]:
    """
    Read models that passed filtering from pairwise_filtering_decisions.csv.
    
    Extracts models by resolving redundancies from the is_distinct column.
    The decisions CSV has format: model1, model2, is_distinct, ...
    
    If pairwise_filtering_decisions.csv doesn't exist, falls back to reading
    McNemarStatsTestsAcrossPRSCalculations_refactored.csv and running filtering.
    """
    import pandas as pd
    from helper.prs_model_filter import filter_redundant_models
    
    # Try reading from pairwise_filtering_decisions.csv first
    decisions_file = os.path.join(scores_path, 'pairwise_filtering_decisions.csv')
    if os.path.exists(decisions_file):
        print(f"      Reading from {os.path.basename(decisions_file)}")
        df = pd.read_csv(decisions_file)
        
        # Extract all unique models
        all_models = set()
        all_models.update(df['model1'].unique())
        all_models.update(df['model2'].unique())
        
        # Build redundancy graph from is_distinct=False pairs
        redundancy_graph = {}
        for _, row in df[df['is_distinct'] == False].iterrows():
            m1, m2 = row['model1'], row['model2']
            redundancy_graph.setdefault(m1, []).append(m2)
            redundancy_graph.setdefault(m2, []).append(m1)
        
        # Resolve redundancies using substring rule (longer name = filtered)
        models_to_keep = set(all_models)
        
        for model, redundant_peers in redundancy_graph.items():
            if model not in models_to_keep:
                continue
            
            for peer in redundant_peers:
                if peer not in models_to_keep:
                    continue
                
                # Substring rule: if one is substring of other, filter the longer one
                if model in peer and model != peer:  # model is substring of peer
                    to_filter = peer
                elif peer in model and model != peer:  # peer is substring of model
                    to_filter = model
                else:  # No substring relationship, use length
                    to_filter = model if len(model) > len(peer) else peer
                
                models_to_keep.discard(to_filter)
                if to_filter == model:
                    break  # This model is filtered, skip remaining peers
        
        models = sorted(models_to_keep)
        
        # Remove 'combined' if present (single-analysis artifact)
        models = [m for m in models if 'combined' not in m.lower()]
        return models
    
    # Fallback: read McNemar CSV and run filtering (creates decisions on the fly)
    mcnemar_file = os.path.join(scores_path, 'McNemarStatsTestsAcrossPRSCalculations_refactored.csv')
    if os.path.exists(mcnemar_file):
        print(f"      No decisions file, filtering from {os.path.basename(mcnemar_file)}")
        mcnemar_df = pd.read_csv(mcnemar_file)
        results = filter_redundant_models(
            stats_results=mcnemar_df,
            model_column='comparison',
            verbose=False
        )
        models = results['models_to_keep']
        # Remove 'combined' if present
        models = [m for m in models if 'combined' not in m.lower()]
        return models
    
    print(f"      WARNING: No filtering files found in {scores_path}")
    return []

def merge_metadata_csvs(
    product_scores: str,
    summed_scores: str,
    combined_scores: str,
) -> tuple[list[str], list[str]]:
    """
    Get the statistically distinct models from separate analysis.
    
    Returns
    -------
    tuple[list[str], list[str]]
        (product_filtered, summed_filtered) - models kept in each analysis
    """
    print("\n[Step 1] Merging FILTERED stats from productEpi/ and summedEpi/...")
    
    # Read which models passed filtering in each analysis
    print("\n  Reading filtered models from separate analyses:")
    print("    productEpi/: ",product_scores)
    product_filtered = _read_models_to_keep_from_filtering(product_scores)
    print(f"      → Kept: {product_filtered}")
    
    print("    summedEpi/: ",summed_scores)
    summed_filtered = _read_models_to_keep_from_filtering(summed_scores)
    print(f"      → Kept: {summed_filtered}")

    return product_filtered, summed_filtered


def stage_prs_score_files(
    product_scores: str,
    summed_scores: str,
    combined_scores: str,
) -> None:
    """
    Build a combined prsScores/ directory so combine_prs and
    calculate_top_features_in_cohort can glob ``*mixed.prs.csv`` files with
    the renamed model identifiers.

    Expected source filename patterns
    ----------------------------------
    Test set  : ``{model}.{n_features}.mixed.prs.csv``
                e.g. epi.200.mixed.prs.csv
    Holdout   : ``{model}.{n_features}.holdout.mixed.prs.csv``
                e.g. epi.200.holdout.mixed.prs.csv

    Both patterns end in ``mixed.prs.csv`` so a single ``*mixed.prs.csv`` glob
    captures them together.  ``{n_features}`` differs between productEpi and
    summedEpi for the same model, so it is **stripped** from the combined
    filename — the model suffix alone prevents collisions:

        epi.200.mixed.prs.csv            → epi_product.mixed.prs.csv
        epi.350.mixed.prs.csv            → epi_summed.mixed.prs.csv
        epi.200.holdout.mixed.prs.csv    → epi_product.holdout.mixed.prs.csv
        epi.350.holdout.mixed.prs.csv    → epi_summed.holdout.mixed.prs.csv
        main.150.mixed.prs.csv           → main.mixed.prs.csv          (shared)
        main.150.holdout.mixed.prs.csv   → main.holdout.mixed.prs.csv  (shared)

    Downstream consumers (combine_prs, calculate_top_features_in_cohort) apply
    their own ``'holdout' not in f`` / ``'holdout' in f`` filters, so test and
    holdout files coexist here without issue.

    Files for shared models (``main``) are copied once from productEpi only.
    """
    combined_prs_scores = os.path.join(combined_scores, "prsScores")
    _make_dirs(combined_prs_scores)

    print("\n[Step 1] Staging combined prsScores/ ...")

    def _copy_with_rename(src_dir: str, suffix: str, skip_shared: bool = False) -> None:
        """
        Copy all ``*mixed.prs.csv`` files from ``src_dir/prsScores/`` into
        ``combined_prs_scores/``, stripping ``{n_features}`` and inserting
        ``suffix`` after the model-name segment.

        Rename rule
        -----------
        Source layout : {model}.{n_features}[.holdout].mixed.prs.csv
          parts[0]    = model      → kept
          parts[1]    = n_features → DROPPED (differs between analyses)
          parts[2:]   = remainder  → kept verbatim

        Result : {model}[{suffix}][.holdout].mixed.prs.csv

        Parameters
        ----------
        src_dir      : Parent directory containing ``prsScores/``.
        suffix       : Model-name suffix, e.g. ``_product`` or ``_summed``.
        skip_shared  : If True, skip files whose stem belongs to SHARED_MODELS.
                       Use when copying the second source to avoid overwriting
                       the already-staged shared file from the first source.
        """
        src_prs_dir = os.path.join(src_dir, "prsScores")
        if not os.path.isdir(src_prs_dir):
            print(f"  WARNING: prsScores/ not found at {src_prs_dir} — skipping.")
            return

        # Glob catches both test and holdout files; downstream scripts filter.
        # Files containing '.FromAll.' are sub-model decompositions
        # (e.g. all.epi+main.FromAll.mixed.prs.csv) and are excluded —
        # only the top-level model files are staged into the combined dir.
        all_files = sorted(glob.glob(os.path.join(src_prs_dir, "*mixed.prs.csv")))
        prs_files = [f for f in all_files if ".FromAll." not in os.path.basename(f)]
        skipped   = len(all_files) - len(prs_files)
        if skipped:
            print(f"  Skipped {skipped} .FromAll. file(s) in {src_prs_dir}")
        for src_file in prs_files:
            basename  = os.path.basename(src_file)

            # Segment layout: {model}.{n_features}[.holdout].mixed.prs.csv
            #   parts[0] = model      e.g. 'epi', 'cardio', 'main'
            #   parts[1] = n_features e.g. '200'  → stripped in combined dir
            #   parts[2:] = remainder e.g. 'mixed.prs.csv' or 'holdout.mixed.prs.csv'
            parts     = basename.split(".")
            stem      = parts[0]
            remainder = ".".join(parts[2:])   # skip parts[1] (n_features)

            if stem in SHARED_MODELS:
                if skip_shared:
                    continue   # already staged from productEpi
                # shared: strip n_features, no suffix
                # main.150.mixed.prs.csv → main.mixed.prs.csv
                new_basename = f"{stem}.{remainder}"
            else:
                # non-shared: strip n_features, add suffix
                # epi.200.mixed.prs.csv → epi_product.mixed.prs.csv
                new_basename = f"{stem}{suffix}.{remainder}"

            dst_file = os.path.join(combined_prs_scores, new_basename)
            if not os.path.exists(dst_file):
                shutil.copy2(src_file, dst_file)
                print(f"  Staged [{suffix.lstrip('_'):>7}]: {new_basename}")
            else:
                print(f"  Already exists (skipped): {new_basename}")

    _copy_with_rename(product_scores, PRODUCT_SUFFIX, skip_shared=False)
    _copy_with_rename(summed_scores,  SUMMED_SUFFIX,  skip_shared=True)


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 – Scale & bin combined PRS
# ─────────────────────────────────────────────────────────────────────────────

def scale_and_bin_combined(
    merged_df: pd.DataFrame,
    prs_list: list[str],
    combined_scores: str,
    combined_figs: str,
) -> tuple[pd.DataFrame, dict]:
    """
    Replicates the scaling / binning pipeline from combine_prs.main() on the
    already-merged combined DataFrame.

    ALL models present in merged_df (including prs_all) are scaled and binned.
    The use_all flag controls only plotting in Step 3, not scaling/binning.

    Returns
    -------
    combined_df     : pd.DataFrame  — merged + scaled + binned training set
    training_stats  : dict          — mean/std per raw PRS column (for holdout)
    """
    print("\n[Step 2] Scaling and binning combined PRS ...")

    df = merged_df.copy()

    # ── Scale all raw PRS columns ─────────────────────────────────────────
    prs_raw_cols = [f"prs_{m}" for m in prs_list if f"prs_{m}" in df.columns]
    scaled_df = scale_data(df[prs_raw_cols])
    scaled_df.columns = [f"scaled_{c}" for c in prs_raw_cols]
    df = df.merge(scaled_df, left_index=True, right_index=True, how="left")

    # ── Bin each PRS model & accumulate prevalence ────────────────────────
    prevalence_df = pd.DataFrame()

    for model in prs_list:
        scaled_col = f"scaled_prs_{model}"
        if scaled_col not in df.columns:
            print(f"  WARNING: {scaled_col} not found — skipping bin step.")
            continue

        prevalence, df_binned = bin_prs(df, scaled_col, model)

        # Visual identity — fall back to defaults if model not in combined map
        prevalence["marker"] = COHORT_MARKERS_COMBINED.get(model, "o")
        prevalence["color"]  = COHORT_COLORS_COMBINED.get(model, "#999999")
        prevalence["model"]  = model

        df = df.merge(df_binned, left_index=True, right_index=True, how="left")
        prevalence_df = pd.concat([prevalence, prevalence_df], ignore_index=True)

    # ── Compute training statistics for holdout scaling ───────────────────
    training_stats = calculate_training_statistics(df)

    # ── Odds-ratio table ──────────────────────────────────────────────────
    odds_ratio_df = calculate_odds_ratio_for_prs(df, "scaled_prs")
    odds_ratio_df.to_csv(
        os.path.join(combined_scores, "combinedORPRSGroups.csv"), index=False
    )

    # ── Prevalence plot ───────────────────────────────────────────────────
    prevalence_df.to_csv(
        os.path.join(combined_scores, "combinedPrevalencePRSGroups.csv"), index=False
    )
    create_optimized_prevalence_plot(prevalence_df, combined_figs, "combined_validation")

    # ── Save main combined file ───────────────────────────────────────────
    combined_prs_path = os.path.join(combined_scores, "combinedPRSGroups.csv")
    df.to_csv(combined_prs_path, index=False)
    print(f"  Written: {combined_prs_path}")
    print(f"  Models scaled/binned: {prs_list}")

    return df, training_stats


def scale_and_bin_holdout(
    product_scores: str,
    summed_scores: str,
    combined_scores: str,
    prs_list: list[str],
    training_stats: dict,
) -> None:
    """
    Merges the holdout datasets from both analyses, applies manual scaling
    using training statistics, and writes the combined holdout file.

    ALL models in prs_list are scaled for the holdout set.
    """
    print("\n[Step 2] Processing combined holdout set ...")

    holdout_product_path = os.path.join(product_scores, "combinedPRSGroups.holdout.csv")
    holdout_summed_path  = os.path.join(summed_scores,  "combinedPRSGroups.holdout.csv")

    if not os.path.exists(holdout_product_path):
        print(f"  WARNING: productEpi holdout not found at {holdout_product_path} — skipping.")
        return
    if not os.path.exists(holdout_summed_path):
        print(f"  WARNING: summedEpi holdout not found at {holdout_summed_path} — skipping.")
        return

    df_h_product = pd.read_csv(holdout_product_path)
    df_h_summed  = pd.read_csv(holdout_summed_path)

    # Drop pre-existing scaled columns; apply the same suffix renaming
    drop_pat = re.compile(r"^(scaled_prs_|bin_)")
    df_h_product = df_h_product[[c for c in df_h_product.columns if not drop_pat.match(c)]]
    df_h_summed  = df_h_summed [[c for c in df_h_summed.columns  if not drop_pat.match(c)]]

    df_h_product.rename(columns=_model_rename_map(PRODUCT_SUFFIX), inplace=True)
    df_h_summed.rename(columns=_model_rename_map(SUMMED_SUFFIX),  inplace=True)

    summed_non_shared = [
        c for c in df_h_summed.columns
        if c.startswith("prs_") and not any(c == f"prs_{m}" for m in SHARED_MODELS)
    ]
    df_holdout = df_h_product.merge(
        df_h_summed[["IID", "PHENOTYPE"] + summed_non_shared],
        on=["IID", "PHENOTYPE"],
        how="inner",
    )

    # Scale holdout using training statistics
    df_holdout_scaled = scale_holdout_data_manually(df_holdout, training_stats)
    out_path = os.path.join(combined_scores, "combinedPRSGroups.holdout.csv")
    df_holdout_scaled.to_csv(out_path, index=False)
    print(f"  Written: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 – filter_statistically_distinct_models
# ─────────────────────────────────────────────────────────────────────────────

def run_filter_distinct_models(
    combined_pheno_data: str,
    use_all: bool,
) -> tuple[list, pd.DataFrame]:
    """
    Identifies statistically distinct models from the combined dataset and
    produces all model-comparison figures.
    
    McNemar/precision/recall statistics are computed for ALL models in
    combinedPRSGroups.csv. The use_all flag controls whether "all v combined"
    comparisons are excluded from the filtering logic.
    
    use_all=False (default): "all v combined" excluded → prs_all won't be kept
    use_all=True: "all v combined" included → prs_all may be kept if distinct
    """
    print("\n[Step 3] Filtering statistically distinct models ...")
    print(f"  use_all={use_all} — controls whether prs_all models can be kept in filtered set")
    models_to_keep, filtered_data = filter_statistically_distinct_models(
        combined_pheno_data,
        use_all=use_all  # passed ONLY for plotting decisions
    )
    print(f"  Distinct models retained: {models_to_keep}")
    return models_to_keep, filtered_data


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 – assign_holdout_risk_bin
# ─────────────────────────────────────────────────────────────────────────────

def run_assign_holdout_bins(
    combined_pheno_data: str,
    use_all: bool,
) -> None:
    """
    Build 1,000-bin risk statistics from ``combinedPRSGroups.filtered.csv``
    (the post-filter test set) and apply them to the holdout set.
    
    Also creates a scatter plot showing prs_combined alongside individual models.

    Reads
    -----
    scores/combinedPRSGroups.filtered.csv   ← training set for bin creation;
                                              prs_columns discovered dynamically
                                              from ``scaled_prs_*`` columns
    scores/combinedPRSGroups.holdout.csv    ← holdout set to be binned

    Writes
    ------
    scores/combinedPRSGroups.holdout.filtered.csv
        Holdout with {model}_centile_bin, {model}_high_risk, combined_centile_bin,
        combined_high_risk, scaled_prs_combined, and prop_case_bins columns added.
    scores/prsScores/prs_statistics_scaled_prs_*_bins.csv
        Per-model bin-statistics CSVs (re-loadable via load_bin_statistics()).
    scores/prsScores/combinedORPRSGroups.holdout.filtered.csv
        Odds-ratio table for the binned holdout PRS.
    figures/prsPlots/
        Prevalence and case-control histogram plots.
    figures/modelComparisons/all_high_risk_cases_with_combined.png
        Scatter plot showing prs_combined with individual models.
    """
    print("\n[Step 4] Assigning holdout risk bins from filtered training set ...")
    combined_scores = os.path.join(combined_pheno_data, "scores")
    combined_figs   = os.path.join(combined_pheno_data, "figures")
    _make_dirs(
        os.path.join(combined_scores, "prsScores"),
        os.path.join(combined_figs,   "prsPlots"),
        os.path.join(combined_figs,   "modelComparisons"),
    )
    run_holdout_bins(
        scoresPath = combined_scores,
        figPath    = combined_figs,
    )
    
    # Create prs_combined scatter plot
    _create_combined_scatter_plot(combined_pheno_data, use_all)


def _create_combined_scatter_plot(
    combined_pheno_data: str,
    use_all: bool,
) -> None:
    """
    Create scatter plot showing prs_combined alongside all filtered models
    from the validation set (HOLDOUT SET).
    
    Uses the plot_holdout_with_combined function that properly
    handles column naming (centile_bin vs bin_) and shows each model with
    correct colors.
    """
    print("\n  Creating HOLDOUT scatter plot with prs_combined ...")
    
    combined_scores = os.path.join(combined_pheno_data, "scores")
    combined_figs = os.path.join(combined_pheno_data, "figures", "modelComparisons")
    
    # Load filtered models from validation set
    filtered_path = os.path.join(combined_scores, "combinedPRSGroups.filtered.csv")
    if not os.path.exists(filtered_path):
        print(f"    WARNING: {filtered_path} not found — skipping plot")
        return
    
    filtered_df = pd.read_csv(filtered_path)
    
    # Get models from filtered validation set
    models_to_keep = [
        c.replace('bin_', '')
        for c in filtered_df.columns
        if c.startswith('bin_')
    ]
    
    # Filter by use_all
    all_variants = {'all', 'all_product', 'all_summed'}
    if not use_all:
        models_to_keep = [m for m in models_to_keep if m not in all_variants]
    
    print(f"    Using {len(models_to_keep)} filtered models from validation set")
    print(f"    Models: {models_to_keep}")
    
    # Load holdout with prs_combined
    holdout_path = os.path.join(combined_scores, "combinedPRSGroups.holdout.filtered.csv")
    if not os.path.exists(holdout_path):
        print(f"    WARNING: {holdout_path} not found — skipping plot")
        return
    
    holdout_df = pd.read_csv(holdout_path)
    
    # Check if prs_combined exists
    if 'scaled_prs_combined' not in holdout_df.columns:
        print("    WARNING: scaled_prs_combined not found — skipping plot")
        return
    
    # Axes are fixed: x=scaled_prs_main, y=scaled_prs_combined
    if 'main' not in models_to_keep:
        print("    WARNING: 'main' not in models_to_keep — skipping holdout scatter")
        return

    try:
        from helper.draw_plots import plot_holdout_with_combined
        import matplotlib.pyplot as plt

        fig, assignments = plot_holdout_with_combined(
            holdout_df=holdout_df,
            models_to_keep=models_to_keep,
            output_path=combined_figs,
            threshold=8,
            figsize=(16, 12)
        )
        
        if fig is not None:
            plt.close(fig)
            
            # Save assignments
            assignments.to_csv(
                os.path.join(combined_scores, 'holdout_case_assignments_with_combined.csv'),
                index=False
            )
            print(f"    ✓ Saved holdout case assignments")
        
    except Exception as e:
        print(f"    ERROR creating holdout scatter plot: {e}")
        import traceback
        traceback.print_exc()
        
# ─────────────────────────────────────────────────────────────────────────────
# Step 5 – Recalculate PRS stats for filtered models
# ─────────────────────────────────────────────────────────────────────────────
      
def recalculate_prs_stats_for_filtered_models(
    combined_pheno_data: str,
    models_to_keep: list[str] | None,
    use_all: bool = False,
) -> None:
    """
    Recalculate PRS performance statistics using only the filtered models from Step 3.
    
    Step 1 merges pre-filtered stats from separate analyses. This step recalculates 
    for the filtered subset from Step 3, which changes comparative metrics:
    
    Unchanged (per-model):
        - TP, FP, TN, FN (confusion matrix)
        - Precision, Recall, Specificity
    
    Changed (comparative):
        - unique_cases: "cases found ONLY by this model" increases when comparison
          set is smaller (fewer competitors)
        - extra_vs_main: may change if filtered set excludes certain models
    
    The original merged stats are backed up to model_recall_precision_improvement.all_models.csv
    
    Parameters
    ----------
    combined_pheno_data : str
        Path to combinedAnalysis directory.
    models_to_keep : list[str] | None
        Filtered models from Step 3. If None, skips recalculation.
    use_all : bool
        Whether to include 'all' model variants.
    """
    if models_to_keep is None:
        print("\n[Step 6] Skipped — no filtered model list available (Step 3 not run)")
        return
  
    print(f"\n[Step 6] Recalculating PRS stats for {len(models_to_keep)} filtered models...")
    print(f"  Models: {', '.join(models_to_keep)}")
  
    combined_scores = os.path.join(combined_pheno_data, 'scores')
    stats_file = os.path.join(combined_scores, 'model_recall_precision_improvement.csv')
  
    # Backup Step 1 stats (all models) if not already done
    backup_file = stats_file.replace('.csv', '.all_models.csv')
    if os.path.exists(stats_file) and not os.path.exists(backup_file):
        import shutil
        shutil.copy2(stats_file, backup_file)
        print(f"    ✓ Backed up Step 1 stats (all models) to {os.path.basename(backup_file)}")
      
    # Recalculate with filtered models
    calculate_prs_stats(
        pheno_data=combined_pheno_data,
        model_type='prs',
        include_all=use_all,
        models_to_keep=models_to_keep,
    )
  
    if os.path.exists(stats_file):
        print(f"    ✓ Recalculated stats saved to {os.path.basename(stats_file)}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 6 – clinical_measure_performance
# ─────────────────────────────────────────────────────────────────────────────

def run_clinical_analysis(
    pheno: str,
    combined_pheno_data: str,
    results_path: str,
) -> None:
    """
    Runs the clinical-measure performance analysis (AUC, NRI, Sankey) on the
    combined holdout PRS set.

    clinicalEnvironmentHoldout.csv does not exist in combinedAnalysis; the
    function resolves it automatically from the sibling productEpi/ or
    summedEpi/ directory (they are identical files).
    """
    print("\n[Step 6] Running clinical measure performance analysis ...")
    # Locate clinicalEnvironmentHoldout.csv in a sibling analysis directory
    import os as _os
    parent = _os.path.dirname(_os.path.abspath(combined_pheno_data))
    clinical_source_path = None
    for candidate in ('productEpi', 'summedEpi'):
        f = _os.path.join(parent, candidate, 'clinicalEnvironmentHoldout.csv')
        if _os.path.exists(f):
            clinical_source_path = f
            print(f'  Using clinicalEnvironmentHoldout.csv from {candidate}/')
            break
    run_clinical_performance(
        pheno=pheno,
        pheno_data=combined_pheno_data,
        results_path=results_path,
        clinical_source_path=clinical_source_path,
    )
  
      
# ─────────────────────────────────────────────────────────────────────────────
# Step 7 – calculate_top_features_in_cohort
# ─────────────────────────────────────────────────────────────────────────────
      
def run_feature_importance(
    combined_pheno_data: str,
    feature_scores_file: str,
    raw_env_features_file : str,
#   threshold: float,
    use_all: bool
) -> None:
    """
    Runs the pipeline for genomic data: SHAP / Random-Forest feature-importance analysis on the combined
    model set, and measured against clinical statistics in pipeline.

    The feature_scores_file should be a concatenation of the filtered feature
    scores from both analyses, with model names updated to include suffixes.
    If the caller has already prepared such a file, it is used directly;
    otherwise a merged version is built on the fly from both analysis
    directories (see _build_combined_feature_scores()).
    raw_features_file is a path to the raw environmental features for indivuals in dataset
    """
    print("\n[Step 5] Calculating feature importance in combined cohorts ...")
    run_pipeline(
          pheno_data=combined_pheno_data,
          feature_scores_file=feature_scores_file,
          raw_features_file=raw_env_features_file,
          use_all=use_all
    )
  
  
def _build_combined_feature_scores(
    product_pheno_data: str,
    summed_pheno_data: str,
    combined_scores: str,
) -> str:
    """
    If a pre-suffixed combined feature-scores file does not yet exist, merge
    the filtered feature-scores CSVs from both analysis directories, applying
    the same model-name suffixing convention.

    Returns the path of the combined feature-scores file.
    """
    out_path = os.path.join(combined_scores, "featureScoresReducedFinalModel.combined.filtered.csv")
  
    if os.path.exists(out_path):
        print(f"  Found existing combined feature scores: {out_path}")
        return out_path
  
    frames = []
    for pheno_data, suffix in [
        (product_pheno_data, PRODUCT_SUFFIX),
        (summed_pheno_data,  SUMMED_SUFFIX),
    ]:
        candidate_paths = [
            os.path.join(pheno_data, "scores", "featureScoresReducedFinalModel.filtered.csv"),
            os.path.join(pheno_data, "scores", "featureScores.filtered.csv"),
        ]
        src = next((p for p in candidate_paths if os.path.exists(p)), None)
        if src is None:
            print(f"  WARNING: No feature-scores file found under {pheno_data}/scores/ — skipping.")
            continue
  
        df = pd.read_csv(src)
        # Suffix non-shared model names in the 'model' column
        if "model" in df.columns:
            df["model"] = df["model"].apply(
                lambda v: _suffix_model_name(str(v), suffix) if pd.notna(v) else v
            )
        frames.append(df)
  
    if not frames:
        raise FileNotFoundError(
            "No feature-scores CSVs found in either analysis directory.  "
            "Provide --feature_scores_file explicitly."
        )
      
    combined = pd.concat(frames, ignore_index=True)
    #keep only one copy of main features
    combined.drop_duplicates(subset=['feature','model'],keep='first',inplace=True)
    combined.to_csv(out_path, index=False)
    print(f"  Written combined feature scores: {out_path}")
    return out_path



# ─────────────────────────────────────────────────────────────────────────────
# Main orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def main(
    pheno_results: str,
    pheno: str,
    results_path: str,
    feature_scores_file: str | None,
    threshold: float,
    steps: set[int],
    raw_env_features_file: str,
    use_all: bool = False,
) -> None:
    """
    Full combined GxG + G+G PRS analysis pipeline.

    Pipeline Steps
    --------------
    1. Merge & stage: 
       a) Merge productEpi + summedEpi data and their PRE-FILTERED stats
       b) Identify statistically distinct models from merged stats
       c) Calculate detailed stats (TP/FP/TN/FN, within-group McNemar) for filtered models
       d) Run cross-analysis McNemar tests (product vs summed) for filtered models
       e) Calculate stats for NEW models (e.g., 'combined' meta-score)
       
       This respects filtering from separate analyses while adding cross-analysis comparisons
       only for models that were statistically distinct in their original pipelines.
       
    2. Scale & bin: Standardize PRS scores, assign bins for training and holdout
    3. Filter models: Identify statistically distinct models using correlation/overlap
                      (further filtering on top of Step 1)
    4. Assign holdout bins: Apply filtered model thresholds to holdout set
    5. Clinical performance: Compare PRS to clinical measures (AUC, NRI)
    6. Feature importance: Identify top features driving each model
    7. Recalculate PRS stats: Update comparative metrics (unique_cases, extra_vs_main)
                              for filtered models from Step 3.

    Parameters
    ----------
    pheno_results       : Root results directory for the phenotype.
                          Expected sub-directories: productEpi/, summedEpi/.
    pheno               : Phenotype label (e.g. 'type2Diabetes').
    results_path        : Path to the broader results directory (for clinical data).
    feature_scores_file : Path to a pre-built combined feature-scores CSV.
                          If None, it is auto-built from both analysis dirs.
    threshold           : Z-score threshold for SHAP feature selection.
    steps               : Which pipeline steps to run (1–7).
    raw_env_features_file   : File path to raw environmental data for all individuals in dataset
    use_all             : Include prs_all models (all_product / all_summed) in
                          scaling, binning, and all downstream steps.
                          Default False.  Changing this flag requires re-running
                          Steps 2–7.
    """

    # ── Resolve paths ─────────────────────────────────────────────────────
    product_pheno = os.path.join(pheno_results, PRODUCT_EPI)
    summed_pheno  = os.path.join(pheno_results, SUMMED_EPI)
    combined_pheno = os.path.join(pheno_results, COMBINED_DIR)

    product_scores = os.path.join(product_pheno,  "scores")
    summed_scores  = os.path.join(summed_pheno,   "scores")
    combined_scores = os.path.join(combined_pheno, "scores")

    combined_figs         = os.path.join(combined_pheno, "figures")
    combined_figs_model   = os.path.join(combined_figs,  "modelComparisons")
    combined_figs_feat    = os.path.join(combined_figs,  "importantCohortFeatures")
    combined_figs_clin    = os.path.join(combined_figs,  "clinicalPerformance")

    _make_dirs(
        combined_scores,
        os.path.join(combined_scores, "importantCohortScores"),
        combined_figs_model,
        combined_figs_feat,
        combined_figs_clin,
    )

    print("=" * 72)
    print(f"  Combined GxG + G+G PRS Analysis — {pheno}")
    print(f"  productEpi  : {product_pheno}")
    print(f"  summedEpi   : {summed_pheno}")
    print(f"  output      : {combined_pheno}")
    print("=" * 72)

    # ── STEP 1 – Merge & stage ─────────────────────────────────────────────
    if 1 in steps:
        # Merge pre-filtered stats from separate analyses
        # Returns which models passed filtering in each analysis
        product_filtered, summed_filtered = merge_metadata_csvs(
            product_scores, summed_scores, combined_scores
        )
      
        merged_df, prs_list = merge_prs_data(
            product_scores, summed_scores, combined_scores,product_filtered,summed_filtered
        )
        stage_prs_score_files(product_scores, summed_scores, combined_scores)
        
        print(f"\n  Merged pre-filtered models:")
        print(f"    productEpi: {product_filtered}")
        print(f"    summedEpi:  {summed_filtered}")
        
        prs_cols  = [c for c in merged_df.columns if c.startswith("prs_") and "scaled" not in c]
        prs_list  = [c.replace("prs_", "") for c in prs_cols]
        
    else:
        # If skipping Step 1, reload from disk for downstream steps
        merged_df = pd.read_csv(os.path.join(combined_scores, "combinedPRSGroups.csv"))
        prs_cols  = [c for c in merged_df.columns if c.startswith("prs_") and "scaled" not in c]
        prs_list  = [c.replace("prs_", "") for c in prs_cols]
        print(f"\n[Step 1] Skipped — using existing combinedPRSGroups.csv ({len(merged_df):,} rows)")

    # ── STEP 2 – Scale & bin ───────────────────────────────────────────────
    if 2 in steps:
      #create combinedPRS files for all models as a starting point
        combined_df, training_stats = scale_and_bin_combined(
            merged_df, prs_list, combined_scores, combined_figs
        )
        scale_and_bin_holdout(
            product_scores, summed_scores, combined_scores,
            prs_list, training_stats
        )
    else:
        combined_df = pd.read_csv(os.path.join(combined_scores, "combinedPRSGroups.csv"))
        training_stats = None
        print(f"\n[Step 2] Skipped — using existing combinedPRSGroups.csv "
              f"({len(combined_df):,} rows)")


    # ── STEP 3 – Filter distinct models ───────────────────────────────────
    if 3 in steps:
        # Calculate stats for ALL models in combined dataset
        print("\n  Calculating stats for all models in combined dataset...")
        
        # Extract models for informational output
        all_models = [c.replace('prs_', '') for c in combined_df.columns 
                      if c.startswith('prs_') and 'scaled' not in c]
      
        if not all_models:
            print("    WARNING: No models found in combined dataset")
        else:
            print(f"    Found {len(all_models)} models in combined dataset: {all_models}")
          
            # Calculate PRS stats for filtered models (TP/FP/TN/FN, McNemar within groups)
            print("\n  Calculating detailed stats for filtered models...")
            calculate_prs_stats(
              pheno_data=combined_pheno,
              model_type='prs',
              include_all=use_all,
              models_to_keep=None,
            )
        
                
            models_to_keep, _ = run_filter_distinct_models(combined_pheno, use_all=use_all)

    else:
        print("\n[Step 3] Skipped — filter_statistically_distinct_models")
        models_to_keep = None

    # ── STEP 4 – Assign holdout risk bins ─────────────────────────────────
    # Requires combinedPRSGroups.filtered.csv produced by Step 3.
    if 4 in steps:
        run_assign_holdout_bins(combined_pheno, use_all=use_all)
    else:
        print("\n[Step 4] Skipped — assign_holdout_risk_bin")

    # ── STEP 5 – Clinical performance ─────────────────────────────────────
    if 5 in steps:
        run_clinical_analysis(pheno, combined_pheno, results_path)
    else:
        print("\n[Step 5] Skipped — clinical_measure_performance")

    # ── STEP 6 – Recalculate PRS stats for filtered models ────────────────
    if 6 in steps:
        # Load models_to_keep if not already loaded from Step 3
        if models_to_keep is None and 3 not in steps:
            filtered_csv = os.path.join(combined_scores, "combinedPRSGroups.filtered.csv")
            if os.path.exists(filtered_csv):
                df_filt = pd.read_csv(filtered_csv, nrows=1)
                models_to_keep = [
                    c.replace('scaled_prs_', '') for c in df_filt.columns
                    if c.startswith('scaled_prs_') 
                    and 'threshold' not in c
                ]
                print(f"  Loaded models_to_keep from filtered CSV: {models_to_keep}")
        recalculate_prs_stats_for_filtered_models(combined_pheno, models_to_keep, use_all)
    else:
        print("\n[Step 6] Skipped — recalculate PRS stats for filtered models")

    print("\n" + "=" * 72)
    print("  Combined pipeline complete.")
    print(f"  All outputs written to: {combined_pheno}")
    print("=" * 72)
    
    # ── STEP 7 – Feature importance ────────────────────────────────────────
    if 7 in steps:
        if feature_scores_file is None:
            feature_scores_file = _build_combined_feature_scores(
                product_pheno, summed_pheno, combined_scores
            )
        run_feature_importance(combined_pheno, feature_scores_file, raw_env_features_file,use_all)
    else:
        print("\n[Step 7] Skipped — calculate_top_features_in_cohort")
# ─────────────────────────────────────────────────────────────────────────────
# CLI entry-point
# ─────────────────────────────────────────────────────────────────────────────

def _parse_steps(raw: str) -> set[int]:
    """Convert '1,2,3' or '1-5' or '3' into a set of integers."""
    steps: set[int] = set()
    for part in raw.split(","):
        part = part.strip()
        if "-" in part:
            start, end = part.split("-")
            steps.update(range(int(start), int(end) + 1))
        else:
            steps.add(int(part))
    return steps


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine GxG (productEpi) and G+G (summedEpi) PRS analyses "
                    "into a unified dataset and run the full downstream pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--pheno_results",
        required=True,
        help="Root results directory containing productEpi/ and summedEpi/ sub-folders.",
    )
    parser.add_argument(
        "--pheno",
        required=True,
        help="Phenotype label (e.g. 'type2Diabetes').",
    )
    parser.add_argument(
        "--results_path",
        required=True,
        help="Broader results directory used by clinical_measure_performance "
             "to locate participant_environment.csv.",
    )
    parser.add_argument(
        "--feature_scores_file",
        default=None,
        help="Path to a pre-built combined feature-scores CSV (with suffixed model names). "
             "If omitted, the script attempts to auto-build one from both analysis directories.",
    )
    parser.add_argument(
        "--use_all",
        action="store_true",
        default=False,
        help="Include prs_all models (all_product/all_summed) in Step 3 plots (scatter plots, correlation matrix). McNemar/precision/recall stats are ALWAYS computed for all models regardless of this flag. Default: False (exclude from plots).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.99,
        help="Z-score threshold for SHAP-based top-feature identification.",
    )
  
    parser.add_argument(
        "--raw_env_features_file",
        default=None,
        help="Path to raw environmental measures used in analysis across all phenotypes.",
    )
    parser.add_argument(
        "--steps",
        default="1,2,3,4,5,6,7",
        help="Comma-separated (or range) list of pipeline steps to run. "
             "Steps: 1=merge+filter+stats+cross_analysis, 2=scale&bin, 3=filter_again, "
             "4=holdout_bins, 5=clinical, 6=features, 7=stats(filtered). "
             "Step 1 merges pre-filtered stats, identifies distinct models, calculates "
             "detailed stats, and runs cross-analysis comparisons. Step 3 does additional "
             "filtering. Step 7 recalculates for Step 3's filtered models.",
    )
  
    args = parser.parse_args()
  
    # Allow environment-variable overrides consistent with the other scripts
    pheno_results      = args.pheno_results or os.environ.get("PHENO_RESULTS")
    pheno              = args.pheno          or os.environ.get("PHENO")
    results_path       = args.results_path   or os.environ.get("RESULTS_PATH")
    feature_scores_file = args.feature_scores_file or os.environ.get("FEATURE_SCORES_FILE_FILTERED")
  
#   pheno_results = '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
#   pheno = 'type2Diabetes'
#   results_path  = '/Users/kerimulterer/prsInteractive/results'
#   feature_scores_file = None
#   use_all = False
#   threshold = 1.99
#   use_all=False
#   steps='1,2,3,4,5,6,7'
  
    if not pheno_results:
        raise ValueError("Provide --pheno_results or set PHENO_RESULTS env var.")
    if not pheno:
        raise ValueError("Provide --pheno or set PHENO env var.")
    if not results_path:
        raise ValueError("Provide --results_path or set RESULTS_PATH env var.")
    if not args.raw_env_features_file:
       raw_env_features_file = f"{results_path}/participant_environment.csv"
      
      
    steps = _parse_steps(args.steps)
#   steps = _parse_steps(steps)
  
  
    print(f"[PYTHON] pheno_results      : {pheno_results}")
    print(f"[PYTHON] pheno              : {pheno}")
    print(f"[PYTHON] results_path       : {results_path}")
    print(f"[PYTHON] feature_scores_file: {feature_scores_file}")
    print(f"[PYTHON] use_all            : {args.use_all}")
    print(f"[PYTHON] threshold          : {args.threshold}")
    print(f"[PYTHON] raw_env_features_file: {args.raw_env_features_file}")
    print(f"[PYTHON] steps              : {sorted(steps)}")
  
    main(
      pheno_results=pheno_results,
      pheno=pheno,
      results_path=results_path,
      feature_scores_file=feature_scores_file,
      raw_env_features_file=raw_env_features_file,
      threshold=args.threshold,
      steps=steps,
      use_all=args.use_all,
    )
  
#   main(
#     pheno_results=pheno_results,
#     pheno=pheno,
#     results_path=results_path,
#     feature_scores_file=feature_scores_file,
#     threshold=threshold,
#     steps=steps,
#     use_all=use_all,
#   )

  