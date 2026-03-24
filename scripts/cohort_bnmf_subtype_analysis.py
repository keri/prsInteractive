#!/usr/bin/env python3
"""
cohort_bnmf_subtype_analysis.py

Performs bNMF subtype analysis for each cohort using individuals and features
derived from the cohort analysis pipeline outputs.

For each cohort the script:
  1. Selects individuals
       High-risk cases   : bin_{cohort} > high_risk_threshold  AND PHENOTYPE == 2
       Low-control ctrls : bin_{cohort} < low_control_threshold AND PHENOTYPE == 1

  2. Selects features from cohort analysis outputs
       Clinical  — FDR-significant features from compare_clinical_features_across_cohorts
                   (values taken from the raw untransformed clinical file)
       Genomic   — SHAP-important features from calculate_top_features_in_cohort
                   (values taken as per-individual weighted contributions from
                    the matching *.mixed.prs.csv model file)

  3. Runs consensus NMF (bNMF) to identify phenotypic subtypes within the
     combined high-risk / low-control population of each cohort.

Inputs
------
  Required:
    pheno_data          — path to phenotype results directory
    raw_features_file   — path to raw (untransformed) clinical features CSV

  Prerequisite outputs from cohort analysis pipeline (run first):
    {pheno_data}/scores/importantCohortScores/
        importantFeaturesAcrossCohortsAndTrainingData.Filtered.{strategy}.csv
    {pheno_data}/scores/cohortClinicalComparison/
        clinical_comparison_all_{use_set}.csv
    {pheno_data}/scores/combinedPRSGroups.holdout.filtered.csv  (or .filtered.csv)
    {pheno_data}/scores/prsScores/*.mixed.prs.csv

Outputs (per cohort under {pheno_data}/scores/cohortBnmf/{cohort}/)
-------
    cluster_assignments.csv   — soft + hard cluster membership per IID
    feature_loadings.csv      — H matrix (feature weights per cluster)
    basis_matrix.csv          — normalised W columns per individual
    cophenetic_curve.csv      — cophenetic correlation per k tested
    cluster_profile.csv       — per-cluster median / IQR statistics

Summary ({pheno_data}/scores/cohortBnmf/)
    cohort_bnmf_overview.csv  — one-row summary per cohort

Figures ({pheno_data}/figures/cohortBnmf/{cohort}/)
    cophenetic_curve.png
    consensus_map_k{k}.png
    cluster_profile.png
    feature_loadings.png

Notes
-----
  Core NMF routines (consensus clustering, cophenetic selection, cluster
  characterisation) are shared with prs_bnmf_subtype_analysis.py via import.
"""

import os
import sys
import glob
import re
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

warnings.filterwarnings("ignore")

# Import shared NMF core from existing bNMF script
import prs_bnmf_subtype_analysis as bnmf_core


# ============================================================================
# CONSTANTS
# ============================================================================

HIGH_RISK_BIN_THRESHOLD  = 8    # bin > this value → high-risk   (deciles 9–10, top 20%)
LOW_CONTROL_BIN_THRESHOLD = 3   # bin < this value → low-control (deciles 1–2, bottom 20%)
                                 # Matches the threshold used in the cohort clinical/genomic
                                 # feature analysis (controls_lt3 selection criterion).
MIN_INDIVIDUALS = 30            # minimum cohort size to attempt bNMF

CLUSTER_PALETTE = bnmf_core.CLUSTER_PALETTE

# Metadata columns present in *.mixed.prs.csv that are not features
PRS_META_COLS = {
    'prs', 'scaled_prs', 'high_risk', 'PHENOTYPE', 'centile_bin',
    'bin', 'combined_high_risk', 'scaled_prs_combined', 'bin_combined',
}


# ============================================================================
# DATA LOADING
# ============================================================================

def load_cohort_genomic_features(scores_path, filter_strategy='tiered'):
    """
    Load SHAP-based genomic features from cohort analysis (step 1 output).

    Returns
    -------
    pd.DataFrame  columns: cohort, feature, matching_z_score, coefs,
                           odds_ratio, [specificity_tier]
                  Empty DataFrame if file not found.
    """
    path = os.path.join(
        scores_path,
        'importantCohortScores',
        f'importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv',
    )
    if not os.path.exists(path):
        print(f"  [WARN] Genomic features file not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    # Normalise column name: feature_set → cohort
    if 'feature_set' in df.columns and 'cohort' not in df.columns:
        df.rename(columns={'feature_set': 'cohort'}, inplace=True)

    n_cohorts = df['cohort'].nunique() if 'cohort' in df.columns else 0
    print(f"  Genomic features loaded : {len(df)} rows, {n_cohorts} cohorts")
    return df


def load_cohort_clinical_features(scores_path, use_set='holdout'):
    """
    Load significant clinical features from cohort analysis (step 2 output).
    Keeps only FDR-significant rows.

    Returns
    -------
    pd.DataFrame  columns: cohort, feature, risk_tier, effect_size_r,
                           p_value_fdr, significant_fdr, ...
                  Empty DataFrame if file not found.
    """
    path = os.path.join(
        scores_path,
        'cohortClinicalComparison',
        f'clinical_comparison_all_{use_set}.csv',
    )
    if not os.path.exists(path):
        print(f"  [WARN] Clinical features file not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    if 'significant_fdr' in df.columns:
        df = df[df['significant_fdr'] == True].copy()

    n_cohorts = df['cohort'].nunique() if 'cohort' in df.columns else 0
    print(f"  Clinical features loaded: {len(df)} significant rows, {n_cohorts} cohorts")
    return df


def load_prs_bins(scores_path, use_set='holdout'):
    """
    Load combined PRS bin assignments for individual selection.

    "Filtered" in the filename means QC-filtered (ancestry, relatedness, etc.),
    NOT risk-filtered.  combinedPRSGroups.holdout.filtered.csv is the standard
    output and contains the full QC-passed holdout population including all
    controls.

    Search order
    ------------
    1. combinedPRSGroups.{holdout.}filtered.csv     — QC-filtered holdout (standard)
    2. combinedPRSGroups.{holdout.}all.csv          — alternative full-pop name
    3. combinedPRSGroups.{holdout.}csv              — no-suffix variant
    4. combinedPRSGroups.filtered.csv / .csv        — non-holdout fallbacks

    Returns
    -------
    pd.DataFrame  indexed by IID, columns include PHENOTYPE, bin_{model}, ...
    """
    holdout_prefix = 'holdout.' if use_set == 'holdout' else ''

    candidates = [
        # Primary: QC-filtered holdout file — standard pipeline output
        os.path.join(scores_path, f'combinedPRSGroups.{holdout_prefix}filtered.csv'),
        # Alternatives with different suffixes
        os.path.join(scores_path, f'combinedPRSGroups.{holdout_prefix}all.csv'),
        os.path.join(scores_path, f'combinedPRSGroups.{holdout_prefix}csv'),
        # Non-holdout fallbacks
        os.path.join(scores_path,  'combinedPRSGroups.filtered.csv'),
        os.path.join(scores_path,  'combinedPRSGroups.all.csv'),
        os.path.join(scores_path,  'combinedPRSGroups.csv'),
    ]
    # Deduplicate while preserving order
    seen = set()
    candidates = [p for p in candidates if not (p in seen or seen.add(p))]

    path = None
    for candidate in candidates:
        if os.path.exists(candidate):
            path = candidate
            break

    if path is None:
        raise FileNotFoundError(
            f"No PRS bins file found in {scores_path}. Tried:\n"
            + "\n".join(f"  {c}" for c in candidates)
        )

    print(f"  PRS bins source         : {os.path.basename(path)}")

    df = pd.read_csv(path)
    if 'IID' not in df.columns and 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'IID'}, inplace=True)
    if 'IID' not in df.columns:
        raise ValueError(
            f"No IID column found in {path}. "
            f"Expected 'IID' or 'Unnamed: 0'. Columns present: {df.columns.tolist()}"
        )
    df = df.set_index('IID')

    # ── Normalize bin columns to a 0-10 scale ─────────────────────────────
    # The pipeline may produce decile (1-10), centile (1-100), or millile
    # (1-1000) bins depending on the scoring script version.  All thresholds
    # in this script assume a 1-10 scale, so we normalise here once.
    bin_cols_all = [c for c in df.columns if c.startswith('bin_')]
    if bin_cols_all:
        # Use the 99th-percentile of observed values (robust to outliers/NaN)
        sample_max = df[bin_cols_all].stack().quantile(0.999)
        if sample_max > 20:          # anything > 20 is not a decile scale
            bin_scale = max(1, round(sample_max / 10))
            df[bin_cols_all] = df[bin_cols_all] / bin_scale
            print(f"  Bin scale detected      : max≈{sample_max:.0f} "
                  f"→ normalised to 0-10 (÷{bin_scale})")
        else:
            print(f"  Bin scale detected      : decile (0-10), no normalisation needed")

    # ── Diagnostic: show PHENOTYPE distribution and controls-in-low-bins ──
    if 'PHENOTYPE' in df.columns:
        ph_counts = df['PHENOTYPE'].value_counts().sort_index()
        print(f"  PRS bins loaded         : {len(df):,} individuals  "
              f"(PHENOTYPE distribution: "
              + ", ".join(f"{k}={v:,}" for k, v in ph_counts.items()) + ")")
        cohort_bin_cols = [c for c in df.columns if c.startswith('bin_') and 'combined' not in c]
        combined_bin_col = next((c for c in df.columns if c.startswith('bin_combined')), None)
        controls = df[df['PHENOTYPE'] == 1]
        if not controls.empty and cohort_bin_cols:
            # Sanity check: how many controls have valid cohort-specific bins?
            for bc in cohort_bin_cols[:5]:   # show first 5 cohorts
                n_valid = controls[bc].notna().sum()
                n_low   = (controls[bc] < 4).sum()
                print(f"    {bc}: {n_valid:,}/{len(controls):,} controls with valid bin, "
                      f"{n_low:,} with bin < 4")
            if combined_bin_col:
                n_valid_comb = controls[combined_bin_col].notna().sum()
                n_low_comb   = (controls[combined_bin_col] < 4).sum()
                print(f"    {combined_bin_col}: {n_valid_comb:,}/{len(controls):,} controls with valid bin, "
                      f"{n_low_comb:,} with bin < 4 (used as fallback for controls with NaN cohort bins)")
        elif controls.empty:
            print("  [WARN] No PHENOTYPE==1 (control) individuals found in PRS bins file.")
    else:
        print(f"  PRS bins loaded         : {len(df):,} individuals  "
              f"(no PHENOTYPE column found)")

    return df


def _extract_model_name(basename):
    """
    Extract model name from a *.mixed.prs.csv filename.

    Examples
    --------
    epi_product.mixed.prs.csv         → epi_product
    main.mixed.prs.holdout.csv        → main
    cardio.mixed.prs.holdout.csv      → cardio
    """
    # Remove trailing .holdout.csv variants first
    name = re.sub(r'\.(holdout\.)?mixed\.prs(\.(holdout))?\.csv$', '',
                  basename, flags=re.IGNORECASE)
    # Fallback: strip any remaining known suffixes
    name = re.sub(r'\.(holdout|mixed|prs|csv)', '', name, flags=re.IGNORECASE)
    return name.strip('._') or basename.split('.')[0]


def load_all_prs_feature_files(scores_path, include_holdout=True, use_set='holdout'):
    """
    Load *.mixed.prs.csv files and return a dict {model_name: DataFrame indexed by IID}.

    File selection is driven by use_set to avoid mixing PRS scores from different
    data partitions (holdout and training individuals have different scaled_prs values):

      use_set='holdout'    → load ONLY files containing 'holdout' in the filename
      use_set='validation' → load ONLY files NOT containing 'holdout' in the filename
      use_set='all'        → load all files; holdout rows stacked on top of training
                             rows, first occurrence per IID kept

    The legacy include_holdout parameter is honoured only when use_set='all'.
    """
    pattern = os.path.join(scores_path, 'prsScores', '*mixed*.csv')
    all_files = glob.glob(pattern)

    if not all_files:
        print(f"  [WARN] No *.mixed.prs.csv files found under {scores_path}/prsScores/")
        return {}

    holdout_files = [f for f in all_files
                     if 'holdout' in os.path.basename(f).lower()]
    train_files   = [f for f in all_files
                     if 'holdout' not in os.path.basename(f).lower()]

    if use_set == 'holdout':
        files_to_use = holdout_files if holdout_files else train_files
        if not holdout_files:
            print("  [WARN] No holdout PRS files found; falling back to training files")
    elif use_set == 'validation':
        files_to_use = train_files
    else:   # 'all'
        files_to_use = train_files + (holdout_files if include_holdout else [])

    model_frames = {}
    for filepath in files_to_use:
        basename = os.path.basename(filepath)
        name = _extract_model_name(basename)
        try:
            df = pd.read_csv(filepath)
            if 'Unnamed: 0' in df.columns:
                df.rename(columns={'Unnamed: 0': 'IID'}, inplace=True)
            if 'IID' not in df.columns:
                print(f"  [WARN] No IID column in {basename}, skipping")
                continue
            df = df.set_index('IID')
            model_frames.setdefault(name, []).append(df)
            print(f"  PRS file loaded         : {basename} → model '{name}'")
        except Exception as exc:
            print(f"  [WARN] Could not load {filepath}: {exc}")

    prs_dict = {}
    for name, frames in model_frames.items():
        if len(frames) == 1:
            prs_dict[name] = frames[0]
        else:
            combined = pd.concat(frames)
            combined = combined[~combined.index.duplicated(keep='first')]
            prs_dict[name] = combined

    print(f"  PRS models available    : {sorted(prs_dict.keys())}")
    return prs_dict


def load_raw_clinical_data(raw_features_file):
    """
    Load raw (untransformed) clinical features.

    Returns
    -------
    pd.DataFrame  indexed by IID, numeric columns only
    """
    df = pd.read_csv(raw_features_file, low_memory=False)
    id_candidates = ['IID', 'Participant ID', 'participant_id', 'eid', 'f.eid']
    id_col = next((c for c in id_candidates if c in df.columns), None)
    if id_col is None:
        raise ValueError(
            f"Cannot find IID column in {raw_features_file}. "
            f"Expected one of: {id_candidates}"
        )
    df = df.rename(columns={id_col: 'IID'}).set_index('IID')
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    print(f"  Raw clinical data loaded: {len(df)} individuals, "
          f"{len(numeric_cols)} numeric features")
    return df[numeric_cols]


# ============================================================================
# INDIVIDUAL SELECTION
# ============================================================================

def _detect_phenotype_coding(prs_bins):
    """
    Detect whether PHENOTYPE is coded as 1/2 (PLINK: control=1, case=2)
    or 0/1 (binary: control=0, case=1).

    Returns
    -------
    int  case_code   — value that means "case"  (2 or 1)
    int  control_code — value that means "control" (1 or 0)
    """
    if 'PHENOTYPE' not in prs_bins.columns:
        return 2, 1   # assume PLINK default

    unique_vals = set(prs_bins['PHENOTYPE'].dropna().unique())
    if 2 in unique_vals:
        # PLINK coding: 1 = control, 2 = case
        return 2, 1
    elif 1 in unique_vals and 0 in unique_vals:
        # Binary coding: 0 = control, 1 = case
        print("  [INFO] PHENOTYPE coded as 0/1 (0=control, 1=case) — adjusting masks.")
        return 1, 0
    else:
        print(f"  [WARN] Unexpected PHENOTYPE values: {sorted(unique_vals)}. "
              f"Defaulting to PLINK coding (1=control, 2=case).")
        return 2, 1


def select_cohort_individuals(
    prs_bins,
    cohort,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    population='both',
):
    """
    Select high-risk cases and/or low-control controls for a given cohort.

    High-risk   : bin_{cohort} > high_risk_threshold  AND PHENOTYPE == case_code
                  (must have a valid cohort-specific bin)
    Low-control : effective_bin < low_control_threshold AND PHENOTYPE == control_code
                  where effective_bin = bin_{cohort} if not NaN, else bin_combined
                  (controls often lack cohort-specific bins because the cohort PRS
                  model may only have been scored for cases; bin_combined is used
                  as a fallback so the full control population is accessible)

    PHENOTYPE coding is detected automatically (PLINK 1/2 or binary 0/1).

    Parameters
    ----------
    population : str  'both' (default) | 'high_risk' | 'low_control'
        Which population sub-group to return.  When 'both', the union of
        high-risk cases and low-control controls is returned (original behaviour).

    Returns
    -------
    pd.Index  IIDs of selected individuals
    dict      {'n_high_risk', 'n_low_control', 'n_total'}
    """
    bin_col = f'bin_{cohort}'
    if bin_col not in prs_bins.columns:
        avail = [c for c in prs_bins.columns if c.startswith('bin_')]
        print(f"  [WARN] '{bin_col}' not found in PRS bins. Available: {avail}")
        return pd.Index([]), {'n_high_risk': 0, 'n_low_control': 0, 'n_total': 0}

    case_code, control_code = _detect_phenotype_coding(prs_bins)

    # ── Bin series for selection ───────────────────────────────────────────────
    # Controls often lack cohort-specific bin assignments (the cohort PRS model
    # may only have been scored for cases / near-threshold individuals). When
    # bin_{cohort} is NaN for a control, fall back to bin_combined if available.
    cohort_bins = prs_bins[bin_col]

    combined_bin_col = next(
        (c for c in prs_bins.columns if c.startswith('bin_combined')), None
    )

    # Build effective bin column: cohort-specific where available, else combined
    if combined_bin_col is not None:
        effective_bins = cohort_bins.combine_first(prs_bins[combined_bin_col])
    else:
        effective_bins = cohort_bins

    # Individuals whose effective bin value is NaN cannot be selected
    valid = effective_bins.notna()

    # High-risk selection always uses the cohort-specific bin (must be valid)
    high_risk_valid = cohort_bins.notna()
    high_risk_mask = (
        high_risk_valid &
        (cohort_bins > high_risk_threshold) &
        (prs_bins['PHENOTYPE'] == case_code)
    )
    # Low-control selection uses the effective bin (cohort → combined fallback)
    low_control_mask = (
        valid &
        (effective_bins < low_control_threshold) &
        (prs_bins['PHENOTYPE'] == control_code)
    )

    n_ctrl_from_cohort_bin   = int((low_control_mask & cohort_bins.notna()).sum())
    n_ctrl_from_combined_bin = int(low_control_mask.sum()) - n_ctrl_from_cohort_bin

    high_risk_iids   = prs_bins[high_risk_mask].index
    low_control_iids = prs_bins[low_control_mask].index

    if population == 'high_risk':
        selected = high_risk_iids
    elif population == 'low_control':
        selected = low_control_iids
    else:  # 'both'
        selected = high_risk_iids.union(low_control_iids)

    stats = {
        'n_high_risk':   len(high_risk_iids),
        'n_low_control': len(low_control_iids),
        'n_total':       len(selected),
        'population':    population,
    }

    if n_ctrl_from_combined_bin > 0:
        print(
            f"    [INFO] {n_ctrl_from_combined_bin:,} low-controls selected via "
            f"bin_combined fallback (cohort-specific bin was NaN for those individuals)"
        )

    # Warn if controls are suspiciously low
    if stats['n_low_control'] == 0:
        n_controls_total = (prs_bins['PHENOTYPE'] == control_code).sum() \
            if 'PHENOTYPE' in prs_bins.columns else 0
        n_bin_valid = valid.sum()
        print(
            f"  [WARN] No low-control individuals found for cohort '{cohort}'.\n"
            f"         Controls in file: {n_controls_total:,}  |  "
            f"  Individuals with valid bin: {n_bin_valid:,}\n"
            f"         Check that PHENOTYPE coding is correct (1=control, 2=case) "
            f"and that bin columns are populated."
        )
    elif stats['n_low_control'] < 50:
        print(
            f"    [WARN] Very few low-controls ({stats['n_low_control']}) for '{cohort}'."
        )

    if population == 'high_risk':
        pop_detail = (f"high-risk (cases): {stats['n_high_risk']:,}  "
                      f"[low-control available: {stats['n_low_control']:,} — not used]")
    elif population == 'low_control':
        pop_detail = (f"low-control (controls): {stats['n_low_control']:,}  "
                      f"[high-risk available: {stats['n_high_risk']:,} — not used]")
    else:
        pop_detail = (f"high-risk: {stats['n_high_risk']:,}  "
                      f"low-control: {stats['n_low_control']:,}")
    print(f"    Individuals — {pop_detail}  →  selected for NMF: {stats['n_total']:,}")
    return selected, stats


# ============================================================================
# FEATURE MATRIX CONSTRUCTION
# ============================================================================

def _find_model_key(cohort, prs_dict):
    """
    Return the key in prs_dict that best matches the cohort name.
    Exact match first, then longest-substring match.
    """
    if cohort in prs_dict:
        return cohort
    # Case-insensitive substring match (longest key wins)
    candidates = sorted(
        prs_dict.keys(),
        key=lambda k: (cohort.lower() in k.lower() or k.lower() in cohort.lower(),
                       len(k)),
        reverse=True,
    )
    for key in candidates:
        if cohort.lower() in key.lower() or key.lower() in cohort.lower():
            return key
    return None


def _filter_prs_dict_to_cohorts(prs_dict, cohorts):
    """
    Keep only PRS model entries whose name EXACTLY matches a valid cohort
    (case-insensitive).

    Substring matching was deliberately removed: it caused composite models
    such as 'epi+main_product' and 'epi+main_summed' to match the 'main'
    cohort (because 'main' is a substring of those names), loading an extra
    two models and triplicating every single SNP from the main cohort.
    """
    cohort_set = {c.lower() for c in cohorts}
    filtered   = {k: v for k, v in prs_dict.items()
                  if k.lower() in cohort_set}
    excluded   = sorted(set(prs_dict.keys()) - set(filtered.keys()))
    if excluded:
        print(f"  Excluded non-cohort PRS models: {excluded}")
    return filtered


def build_cohort_feature_matrix(
    cohort,
    selected_iids,
    genomic_features_df,
    clinical_features_df,
    prs_dict,
    raw_clinical_df,
    min_effect_size=0.0,
    specificity_tier_max=3,
    include_raw_clinical=True,
):
    """
    Build the feature matrix for bNMF for a single cohort.

    Combines:
      Clinical features — significant for this cohort (raw values from
                          raw_clinical_df), prefixed 'clin_'
      Genomic features  — ALL per-variant weighted contribution columns from
                          the cohort's PRS model (prs_dict[model]), prefixed
                          'gen_'.  Every SNP-level column is included; near-
                          constant and highly sparse features are removed by
                          the variance/sparsity filter in run_cohort_bnmf().

    Parameters
    ----------
    cohort                : str
    selected_iids         : pd.Index of IIDs to include
    genomic_features_df   : output of load_cohort_genomic_features()
                            (retained for API compatibility; not used for
                            genomic feature selection in per-cohort mode)
    clinical_features_df  : output of load_cohort_clinical_features()
    prs_dict              : output of load_all_prs_feature_files()
    raw_clinical_df       : output of load_raw_clinical_data()
    min_effect_size       : float, minimum |effect_size_r| filter for clinical
    specificity_tier_max  : int, not used in per-cohort mode (kept for API
                            compatibility with combined mode)

    Returns
    -------
    pd.DataFrame  individuals × features (unscaled, may contain NaN)
    dict          {'n_clinical', 'n_genomic', 'n_individuals'}
    """
    frames = []
    feat_counts = {'n_clinical': 0, 'n_genomic': 0}

    # ------------------------------------------------------------------
    # 1. Clinical features
    # ------------------------------------------------------------------
    if not clinical_features_df.empty and 'cohort' in clinical_features_df.columns:
        cohort_clin = clinical_features_df[
            clinical_features_df['cohort'] == cohort
        ].copy()

        if min_effect_size > 0 and 'effect_size_r' in cohort_clin.columns:
            cohort_clin = cohort_clin[
                cohort_clin['effect_size_r'].abs() >= min_effect_size
            ]

        clin_features = cohort_clin['feature'].unique().tolist()
        available = [f for f in clin_features if f in raw_clinical_df.columns]

        if available:
            clin_sub = raw_clinical_df.loc[
                raw_clinical_df.index.isin(selected_iids), available
            ].copy()
            clin_sub.columns = [f'clin_{c}' for c in clin_sub.columns]
            frames.append(clin_sub)
            feat_counts['n_clinical'] = len(available)
            n_missing = len(clin_features) - len(available)
            print(f"    Clinical  : {len(available)} features used"
                  + (f" ({n_missing} not in data)" if n_missing else ""))
        else:
            print(f"    [WARN] No clinical features found in data for cohort '{cohort}'")

    # ------------------------------------------------------------------
    # 2. Genomic features — ALL per-variant columns from the PRS model
    #    Uses every SNP-level weighted contribution column (excludes meta
    #    columns and principal components).  No SHAP pre-selection filter;
    #    the variance/sparsity filter inside run_cohort_bnmf() removes
    #    near-constant / highly-sparse features before NMF.
    # ------------------------------------------------------------------
    model_key = _find_model_key(cohort, prs_dict)
    if model_key is None:
        print(f"    [WARN] No PRS model matched cohort '{cohort}'. "
              f"Available models: {list(prs_dict.keys())}")
    else:
        prs_df = prs_dict[model_key]
        prs_feature_cols = [
            c for c in prs_df.columns
            if c not in PRS_META_COLS and not c.startswith('PC')
        ]

        if prs_feature_cols:
            gen_sub = prs_df.loc[
                prs_df.index.isin(selected_iids), prs_feature_cols
            ].copy()
            gen_sub.columns = [f'gen_{c}' for c in gen_sub.columns]
            frames.append(gen_sub)
            feat_counts['n_genomic'] = len(prs_feature_cols)
            print(f"    Genomic   : {len(prs_feature_cols)} per-variant features "
                  f"(all columns, model '{model_key}')")
        else:
            print(f"    [WARN] No per-variant feature columns found for "
                  f"cohort '{cohort}' in model '{model_key}'")

    # ------------------------------------------------------------------
    # 3. Raw clinical supplement / fallback
    #    Add all raw clinical features not already covered above.
    #    Clinical features have real values for essentially all individuals
    #    and break the rank-1 sparsity pattern of the genomic-only matrix.
    # ------------------------------------------------------------------
    if include_raw_clinical and raw_clinical_df is not None and not raw_clinical_df.empty:
        already_clin = {c[len('clin_'):]
                        for frm in frames
                        for c in frm.columns if c.startswith('clin_')}
        extra_cols = [c for c in raw_clinical_df.columns if c not in already_clin]
        if extra_cols:
            extra_block = raw_clinical_df.reindex(selected_iids)[extra_cols].copy()
            extra_block.columns = [f'clin_{c}' for c in extra_cols]
            frames.append(extra_block)
            tag = ('fallback (no comparison file)'
                   if feat_counts['n_clinical'] == 0 else 'supplement')
            n_covered = extra_block.notna().any(axis=1).sum()
            print(f"    Clinical (raw {tag}): {len(extra_cols)} features, "
                  f"{n_covered}/{len(extra_block)} individuals have ≥1 value")
            feat_counts['n_clinical'] += len(extra_cols)

    # ------------------------------------------------------------------
    # 4. Merge all feature blocks
    # ------------------------------------------------------------------
    if not frames:
        print(f"    [ERROR] No features available for cohort '{cohort}'")
        return pd.DataFrame(), feat_counts

    merged = frames[0]
    for f in frames[1:]:
        merged = merged.join(f, how='outer')

    # Restrict to selected individuals
    merged = merged.loc[merged.index.isin(selected_iids)]
    feat_counts['n_individuals'] = len(merged)

    print(f"    Matrix    : {merged.shape[0]:,} × {merged.shape[1]} "
          f"({feat_counts['n_clinical']} clinical, {feat_counts['n_genomic']} genomic)")
    return merged, feat_counts


# ============================================================================
# VISUALISATION
# ============================================================================

def _save_fig(fig, path):
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved figure → {os.path.basename(path)}")


def plot_cophenetic_curve(coph_df, output_path, cohort):
    """Line plot of cophenetic correlation vs k.

    Non-collapsed k values are plotted normally.  Collapsed k values (NaN
    cophenetic) are shown as red X markers at y=0 so the user can see which
    k values failed rather than having them silently absent from the plot.
    """
    if coph_df.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4))

    valid   = coph_df.dropna(subset=['cophenetic_correlation'])
    invalid = coph_df[coph_df['cophenetic_correlation'].isna()]

    if not valid.empty:
        ax.plot(valid['k'], valid['cophenetic_correlation'],
                marker='o', color='#1f77b4', linewidth=2, markersize=6,
                label='Stable')

    if not invalid.empty:
        ax.scatter(invalid['k'], [0] * len(invalid),
                   marker='x', color='red', s=80, zorder=5,
                   label='Collapsed (NaN)')
        ax.legend(fontsize=8)

    # Annotate collapse fraction if available
    if 'mean_collapse_frac' in coph_df.columns:
        for _, row in coph_df.iterrows():
            y = row['cophenetic_correlation'] if not np.isnan(row['cophenetic_correlation']) else 0
            ax.annotate(f"{row['mean_collapse_frac']:.0%}",
                        xy=(row['k'], y), xytext=(0, 6),
                        textcoords='offset points', ha='center', fontsize=7,
                        color='red' if np.isnan(row['cophenetic_correlation']) else 'grey')

    ax.set_xlabel('Number of clusters (k)')
    ax.set_ylabel('Cophenetic correlation')
    ax.set_title(f'Consensus NMF Stability — {cohort}')
    ax.set_xticks(coph_df['k'])
    ax.set_ylim(bottom=-0.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    _save_fig(fig, output_path)


def plot_consensus_map(C, k, output_path, cohort):
    """Heatmap of the soft consensus co-clustering matrix."""
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(C, cmap='viridis', ax=ax,
                xticklabels=False, yticklabels=False,
                cbar_kws={'label': 'Co-clustering frequency'})
    ax.set_title(f'Consensus Matrix  k={k} — {cohort}')
    plt.tight_layout()
    _save_fig(fig, output_path)


def plot_cluster_profile(profile_df, output_path, cohort, top_n=20):
    """
    Heatmap of per-cluster feature profiles.

    After QuantileTransformer normalization all raw medians are ~0.5, so
    raw medians are uninformative.  Instead we z-score each feature's
    cluster medians across clusters (z = (median_k - mean_k) / std_k).
    This shows which features are *relatively* higher or lower in each
    cluster compared with the overall average, regardless of the absolute
    scaled value.

    Feature selection uses between-cluster variance of raw medians (computed
    before z-scoring) so we pick the features that actually differ between
    clusters.  Annotations show the z-score and raw median in brackets so
    the reader can see both the direction and the absolute level.
    """
    if profile_df.empty:
        return
    median_cols = [c for c in profile_df.columns
                   if c.endswith('_median') and c not in ('n_median',)]
    if not median_cols:
        return

    # Build raw-median matrix (clusters × features)
    if profile_df.index.name == 'cluster':
        profile_matrix = profile_df[median_cols].copy()
    elif 'cluster' in profile_df.columns:
        profile_matrix = profile_df.set_index('cluster')[median_cols].copy()
    else:
        profile_matrix = profile_df[median_cols].copy()
    profile_matrix.columns = [c.replace('_median', '') for c in median_cols]

    # Select top features by between-cluster variance of raw medians
    if len(profile_matrix.columns) > top_n:
        top_feats = profile_matrix.var(axis=0).nlargest(top_n).index.tolist()
        profile_matrix = profile_matrix[top_feats]

    # Z-score each feature across clusters so all features are on the same
    # scale and the plot shows relative enrichment / depletion per cluster.
    feat_mean = profile_matrix.mean(axis=0)
    feat_std  = profile_matrix.std(axis=0).replace(0, np.nan)
    z_matrix  = (profile_matrix - feat_mean) / feat_std
    z_matrix  = z_matrix.fillna(0)   # constant features → z=0 (not discriminating)

    # Annotation: "z\n(raw)" per cell  e.g.  "+1.23\n(0.61)"
    annot = pd.DataFrame(index=z_matrix.index, columns=z_matrix.columns, dtype=object)
    for feat in z_matrix.columns:
        for clust in z_matrix.index:
            z   = z_matrix.loc[clust, feat]
            raw = profile_matrix.loc[clust, feat]
            annot.loc[clust, feat] = f'{z:+.2f}\n({raw:.2f})'

    # Clip z for colormap symmetry (outliers can dominate otherwise)
    z_clipped = z_matrix.clip(-3, 3)

    fig_h = max(4, len(profile_matrix.columns) * 0.45)
    fig_w = max(4, len(profile_matrix) * 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        z_clipped.T,
        cmap='RdBu_r', center=0, vmin=-3, vmax=3,
        ax=ax,
        annot=annot.T, fmt='',
        annot_kws={'size': 6},
        linewidths=0.5,
        cbar_kws={'label': 'Z-score of cluster median (raw median in brackets)'},
    )
    ax.set_title(
        f'Cluster Profiles — {cohort}\n'
        f'(top {len(profile_matrix.columns)} features by between-cluster variance; '
        f'colour = z-score across clusters)',
        fontsize=9,
    )
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Feature')
    plt.tight_layout()
    _save_fig(fig, output_path)


def plot_feature_loadings(H_df, output_path, cohort, top_n=20):
    """
    Per-cluster horizontal bar charts of cluster-specific feature loadings.

    Raw NMF H-matrix values measure absolute feature importance per cluster
    but features shared across all clusters look equally uninformative even
    if they are the dominant signal.  We z-score each feature's loading
    across clusters (z = (H_k - mean_H) / std_H) so bars show how much
    more (positive z, red) or less (negative z, blue) each feature loads in
    this cluster relative to the cross-cluster average.

    Top N features per cluster are selected by absolute z-score so features
    that truly distinguish this cluster rise to the top.  Clinical features
    are outlined with a darker border to distinguish them from genomic features
    that share the same hue.
    """
    if H_df.empty:
        return

    # Z-score each feature across clusters
    feat_mean = H_df.mean(axis=0)
    feat_std  = H_df.std(axis=0).replace(0, np.nan)
    H_z = ((H_df - feat_mean) / feat_std).fillna(0)

    n_clusters = len(H_df)
    fig, axes = plt.subplots(
        1, n_clusters,
        figsize=(n_clusters * 5, max(4, min(top_n * 0.38, 14))),
    )
    if n_clusters == 1:
        axes = [axes]

    for ax, cluster_label in zip(axes, H_z.index):
        z_row = H_z.loc[cluster_label]
        # Top N by absolute z-score — features that most distinguish this cluster
        top = z_row.abs().nlargest(top_n)
        top_z = z_row[top.index]   # preserve sign

        colours = []
        for f in top_z.index:
            is_clin = f.startswith('clin_')
            if top_z[f] >= 0:
                colours.append('#d62728' if is_clin else '#1f77b4')   # red / blue (enriched)
            else:
                colours.append('#f4a3a3' if is_clin else '#aec7e8')   # pale (depleted)

        labels = [f.replace('clin_', 'C: ').replace('gen_', 'G: ')
                  for f in top_z.index]

        bars = ax.barh(range(len(top_z)), top_z.values, color=colours, alpha=0.88)
        # Darker edge for clinical features
        for bar, f in zip(bars, top_z.index):
            if f.startswith('clin_'):
                bar.set_edgecolor('#7f0000')
                bar.set_linewidth(0.8)
        ax.axvline(0, color='black', linewidth=0.6, linestyle='--')
        ax.set_yticks(range(len(top_z)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_title(cluster_label, fontsize=9)
        ax.set_xlabel('Z-score of NMF loading')

    legend_handles = [
        mpatches.Patch(facecolor='#d62728', edgecolor='#7f0000', linewidth=0.8,
                       label='Clinical (enriched)'),
        mpatches.Patch(facecolor='#f4a3a3', edgecolor='#7f0000', linewidth=0.8,
                       label='Clinical (depleted)'),
        mpatches.Patch(facecolor='#1f77b4', label='Genomic (enriched)'),
        mpatches.Patch(facecolor='#aec7e8', label='Genomic (depleted)'),
    ]
    axes[-1].legend(handles=legend_handles, loc='lower right', fontsize=7)
    plt.suptitle(
        f'Cluster-specific Feature Loadings (z-scored across clusters) — {cohort}',
        fontsize=10,
    )
    plt.tight_layout()
    _save_fig(fig, output_path)


# ============================================================================
# COMBINED bNMF — HELPERS
# ============================================================================

def build_cohort_membership_flags(
    prs_bins,
    cohorts,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
):
    """
    Build a boolean flag matrix for cohort membership.

    Returns
    -------
    pd.DataFrame  indexed by IID
        Columns: is_high_risk_{cohort}, is_low_control_{cohort}  for each cohort.
    """
    flag_frames = []
    for cohort in cohorts:
        col = f'bin_{cohort}'
        if col not in prs_bins.columns:
            continue
        is_case    = prs_bins['PHENOTYPE'] == 2
        is_control = prs_bins['PHENOTYPE'] == 1
        flags = pd.DataFrame({
            f'is_high_risk_{cohort}':   (prs_bins[col] > high_risk_threshold)  & is_case,
            f'is_low_control_{cohort}': (prs_bins[col] < low_control_threshold) & is_control,
        }, index=prs_bins.index)
        flag_frames.append(flags)

    if not flag_frames:
        return pd.DataFrame(index=prs_bins.index)

    return pd.concat(flag_frames, axis=1).fillna(False)


def build_deduplicated_feature_matrix(
    cohorts,
    prs_bins,
    genomic_features_df,
    clinical_features_df,
    prs_dict,
    raw_clinical_df,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    min_effect_size=0.0,
    specificity_tier_max=3,
    include_prs_with_genomic=False,
    population='both',
    include_raw_clinical=True,
):
    """
    Build a combined feature matrix across all cohorts without redundant columns.

    Column structure
    ----------------
    {cohort}_prs    — scaled_prs from each cohort's model.  Added ONLY for cohorts
                      that have NO selected genomic (gen_) features, so the NMF still
                      has a signal for those cohorts.  For cohorts that DO have gen_
                      features, the PRS anchor is suppressed by default because:
                        PRS ≈ Σ(weighted SNP contributions)
                      Including both the aggregate (PRS) and its components (gen_*)
                      introduces near-perfect multicollinearity that NMF absorbs into
                      a spurious dominant component (overfitting / cluster collapse).
                      Set include_prs_with_genomic=True to override this behaviour.

    clin_{feature}  — raw clinical value; ONE column per feature regardless of how many
                      cohorts flag it; weight = max(|effect_size_r|) across cohorts.

    gen_{feature}   — per-individual weighted SNP contribution from the cohort model
                      with the highest matching_z_score for that feature;
                      weight = max z-score.  Replaces the PRS anchor for cohorts
                      that have selected features (see above).

    Individuals are determined by the population parameter:
      'both'        — union of high-risk cases and low-control controls (default)
      'high_risk'   — only high-risk cases
      'low_control' — only low-risk controls
      'separate'    — same as 'both' for the combined matrix (run two separate
                      combined analyses via run_combined_bnmf if needed)

    Returns
    -------
    combined_df      : pd.DataFrame  (union_IIDs × deduplicated features)
    feature_metadata : dict  {column_name: {'cohort_origin', 'type', 'weight'}}
    per_cohort_iids  : dict  {cohort: pd.Index of selected IIDs}
    """
    # 'separate' has no meaning for a single combined matrix; treat as 'both'
    _pop = 'both' if population == 'separate' else population

    per_cohort_iids = {}
    all_iids        = pd.Index([])

    # Collect per-cohort individual sets and union of all IIDs
    for cohort in cohorts:
        iids, _ = select_cohort_individuals(
            prs_bins, cohort,
            population=_pop,
            high_risk_threshold=high_risk_threshold,
            low_control_threshold=low_control_threshold,
        )
        per_cohort_iids[cohort] = iids
        all_iids = all_iids.union(iids)

    if all_iids.empty:
        return pd.DataFrame(), {}, per_cohort_iids

    feature_metadata = {}
    frames           = []

    # ── 1. Genomic features — ALL per-variant columns from ALL models ──────
    #
    # Column deduplication rules (applied BEFORE combine_first):
    #
    #   a) HLA loci  — raw feature name contains '_' (e.g. DRB1_04, A_101).
    #      The same allele may be selected by multiple models with slightly
    #      different coefficients.  Biologically it is the SAME feature, so
    #      we deduplicate: all models use the same gen_{feat} column name and
    #      combine_first() fills real values from whichever model has them,
    #      prioritising the model with the most features.
    #
    #   b) Non-HLA SNPs that appear in ONLY ONE model — no ambiguity, kept
    #      as gen_{feat}.
    #
    #   c) SNP pairs (two rs-IDs joined by '_') that appear in 2+ models —
    #      same pair but different coefficients across product vs summed models.
    #      These are renamed gen_{feat}_{cohort_suffix} using the full cohort
    #      name, e.g. gen_rs123_rs456_epi_prod / gen_rs123_rs456_cardio_sum.
    #      Single SNPs (no '_') should only appear in the 'main' model after
    #      composite models (epi+main_*) are excluded by the cohort filter;
    #      they receive no suffix → gen_rs826962.
    #
    # Because *.mixed.prs.csv files contain ALL genotyped individuals (not
    # just those high-risk in that model), every individual in all_iids has a
    # real value in each model file — NaN→0 imputation is needed only for the
    # rare case where an individual is absent from a particular file entirely.
    # ──────────────────────────────────────────────────────────────────────

    def _model_suffix(model_key: str) -> str:
        """Return the full cohort-qualified suffix, e.g. _epi_prod, _cardio_sum, _main.

        Using the full cohort name avoids ambiguity: _prod alone could mean
        epi_product or cardio_product.  We compress only the verbose trailing
        words 'product'→'prod' and 'summed'→'sum' so the result is short but
        still unambiguous (e.g. epi_product → _epi_prod, cardio_summed → _cardio_sum).
        """
        mk = model_key.lower()
        mk = mk.replace('_product', '_prod').replace('_summed', '_sum')
        return f'_{mk}'   # e.g. _epi_prod, _cardio_sum, _main

    # Build a lookup: feature_name → best z_score (for optional weighting)
    feat_z_scores: dict[str, float] = {}
    if not genomic_features_df.empty and 'cohort' in genomic_features_df.columns:
        z_col = 'matching_z_score' \
                if 'matching_z_score' in genomic_features_df.columns else None
        if z_col:
            for _, row in genomic_features_df.iterrows():
                f = row.get('feature')
                z = abs(float(row[z_col])) if pd.notna(row.get(z_col)) else 0.0
                if f and z > feat_z_scores.get(f, 0.0):
                    feat_z_scores[f] = z

    # Sort models: most features first so combine_first prioritises the
    # richest model's real values for shared HLA columns.
    sorted_models = sorted(
        prs_dict.items(),
        key=lambda kv: -len([c for c in kv[1].columns
                              if c not in PRS_META_COLS and not c.startswith('PC')])
    )

    # Pre-scan: how many models contain each raw feature name?
    from collections import Counter
    feat_model_count: Counter = Counter()
    model_feat_sets: dict[str, set] = {}
    for model_key, prs_df in sorted_models:
        feat_set = {c for c in prs_df.columns
                    if c not in PRS_META_COLS and not c.startswith('PC')}
        model_feat_sets[model_key] = feat_set
        feat_model_count.update(feat_set)

    # Features appearing in 2+ models AND lacking '_' → need model suffix
    snp_needs_suffix: set[str] = {
        f for f, cnt in feat_model_count.items()
        if cnt > 1 and '_' not in f
    }
    if snp_needs_suffix:
        print(f"    Disambiguation: {len(snp_needs_suffix)} non-HLA SNP(s) appear in "
              f"multiple models → annotated with _prod/_sum/_<model> suffix")

    gen_frames: list[tuple[str, pd.DataFrame]] = []
    cohorts_with_genomic: set[str] = set()

    for model_key, prs_df in sorted_models:
        feature_cols = list(model_feat_sets[model_key])
        if not feature_cols:
            continue

        suffix = _model_suffix(model_key)
        blk = prs_df.loc[prs_df.index.isin(all_iids), feature_cols].copy()

        # Rename each column according to the rules above
        new_names = []
        for c in blk.columns:
            if c in snp_needs_suffix:
                # Non-HLA duplicate → add model-type suffix
                new_names.append(f'gen_{c}{suffix}')
            else:
                # HLA locus or unique SNP → no suffix; HLA deduped by combine_first
                new_names.append(f'gen_{c}')
        blk.columns = new_names
        gen_frames.append((model_key, blk))

        # Track cohorts covered by gen_ features (suppresses redundant PRS anchor)
        for cohort in cohorts:
            if _find_model_key(cohort, {model_key: prs_df}) == model_key:
                cohorts_with_genomic.add(cohort)

    if gen_frames:
        # combine_first across models:
        #   • HLA columns (same name across models) → fills real values from
        #     richer model first; no duplication.
        #   • Non-HLA disambiguated columns (gen_feat_prod, gen_feat_sum) →
        #     different names, so they accumulate as separate columns.
        #   • Unique SNP columns → simply appended.
        gen_block = gen_frames[0][1]
        for _, blk in gen_frames[1:]:
            gen_block = gen_block.combine_first(blk)

        frames.append(gen_block)

        n_models   = len(gen_frames)
        n_feats    = gen_block.shape[1]
        n_complete = int(gen_block.notna().all(axis=1).sum())
        print(f"    Genomic   : {n_feats} per-variant features across {n_models} model(s) "
              f"(HLA loci deduplicated; non-HLA duplicates suffixed)")
        print(f"               {n_complete}/{len(all_iids)} individuals have real values "
              f"in every feature; remainder filled by combine_first across models.")

        for col in gen_block.columns:
            # Strip gen_ prefix and any trailing _prod/_sum/_<model> suffix
            # to look up the base feature name in feat_z_scores.
            raw_feat = col[4:]   # remove 'gen_'
            base     = raw_feat  # default: no suffix to strip
            for mk, _ in sorted_models:
                sfx = _model_suffix(mk)
                if raw_feat.endswith(sfx):
                    base = raw_feat[: -len(sfx)]
                    break
            feature_metadata[col] = {
                'cohort_origin': 'combined',
                'type':          'genomic',
                'weight':        feat_z_scores.get(base, 1.0),
            }

        if cohorts_with_genomic and not include_prs_with_genomic:
            print(f"    PRS anchor suppressed for cohorts with gen_ features "
                  f"({sorted(cohorts_with_genomic)}) — gen_ features replace "
                  f"PRS to avoid sum+components redundancy. "
                  f"Pass include_prs_with_genomic=True to override.")

    # ── 2. PRS anchors — added ONLY for cohorts not covered by gen_ features.
    #       Each cohort's scaled_prs is a genuine non-redundant signal for
    #       cohorts where no important SNPs were individually selected.
    # ──────────────────────────────────────────────────────────────────────
    for cohort in cohorts:
        if not include_prs_with_genomic and cohort in cohorts_with_genomic:
            # gen_ features already carry this cohort's genomic signal
            continue
        model_key = _find_model_key(cohort, prs_dict)
        if model_key is None:
            continue
        prs_df  = prs_dict[model_key]
        prs_col = 'scaled_prs' if 'scaled_prs' in prs_df.columns else \
                  ('prs' if 'prs' in prs_df.columns else None)
        if prs_col is None:
            continue
        col_name   = f'{cohort}_prs'
        prs_series = prs_df.loc[prs_df.index.isin(all_iids), prs_col] \
                            .rename(col_name)
        frames.append(prs_series.to_frame())
        feature_metadata[col_name] = {
            'cohort_origin': cohort,
            'type':          'prs',
            'weight':        1.0,
        }

    prs_anchors_added = [c for c, m in feature_metadata.items() if m['type'] == 'prs']
    if prs_anchors_added:
        print(f"    PRS anchors: {len(prs_anchors_added)} cohorts without gen_ features "
              f"({[c.replace('_prs','') for c in prs_anchors_added]})")

    # ── 3. Clinical features — deduplicated, weight = max effect size ──────
    # Accumulate best (max |effect_size_r|) across all cohorts per feature.
    clin_best: dict[str, dict] = {}   # feature_name → {effect_size_r, cohort}
    if not clinical_features_df.empty and 'cohort' in clinical_features_df.columns:
        for cohort in cohorts:
            c_sub = clinical_features_df[
                clinical_features_df['cohort'] == cohort
            ].copy()
            if min_effect_size > 0 and 'effect_size_r' in c_sub.columns:
                c_sub = c_sub[c_sub['effect_size_r'].abs() >= min_effect_size]
            eff_col = 'effect_size_r' if 'effect_size_r' in c_sub.columns else None
            for f in c_sub['feature'].unique():
                if f not in raw_clinical_df.columns:
                    continue
                eff = abs(float(c_sub.loc[c_sub['feature'] == f, eff_col].iloc[0])) \
                    if eff_col else 0.0
                if f not in clin_best or eff > clin_best[f]['weight']:
                    clin_best[f] = {'weight': eff, 'cohort_origin': cohort}

    if clin_best:
        clin_block = raw_clinical_df.loc[
            raw_clinical_df.index.isin(all_iids), list(clin_best.keys())
        ].copy()
        clin_block.columns = [f'clin_{f}' for f in clin_block.columns]
        frames.append(clin_block)
        for f, meta in clin_best.items():
            feature_metadata[f'clin_{f}'] = {
                'cohort_origin': meta['cohort_origin'],
                'type':          'clinical',
                'weight':        meta['weight'],
            }
        print(f"    Clinical  : {len(clin_best)} deduplicated features "
              f"(max |effect_size_r| across cohorts)")

    # ── Raw clinical supplement / fallback ─────────────────────────────────
    # Add ALL raw clinical features not already covered by the comparison file.
    # Clinical features have real values for essentially all individuals —
    # no combine_first imputation artifacts — so they provide genuine biological
    # signal for cluster separation even when the genomic block is sparse.
    # This is the primary fallback when clinical_comparison file is missing.
    if include_raw_clinical and raw_clinical_df is not None and not raw_clinical_df.empty:
        already_clin = {f[len('clin_'):] for f in feature_metadata if f.startswith('clin_')}
        extra_cols = [c for c in raw_clinical_df.columns if c not in already_clin]
        if extra_cols:
            extra_block = raw_clinical_df.reindex(all_iids)[extra_cols].copy()
            extra_block.columns = [f'clin_{c}' for c in extra_cols]
            frames.append(extra_block)
            for c in extra_cols:
                feature_metadata[f'clin_{c}'] = {
                    'cohort_origin': 'all',
                    'type':          'clinical',
                    'weight':        1.0,
                }
            n_covered = extra_block.notna().any(axis=1).sum()
            tag = 'fallback (no comparison file)' if not clin_best else 'supplement'
            print(f"    Clinical (raw {tag}): {len(extra_cols)} features, "
                  f"{n_covered}/{len(extra_block)} individuals have ≥1 value")

    if not frames:
        return pd.DataFrame(), {}, per_cohort_iids

    combined_df = frames[0]
    for frm in frames[1:]:
        combined_df = combined_df.join(frm, how='outer')

    # Re-index to ensure all union IIDs are present (some may only appear in
    # the PRS anchor blocks, not in clinical/genomic slices)
    combined_df = combined_df.reindex(all_iids)

    n_prs  = sum(1 for k, v in feature_metadata.items() if v['type'] == 'prs')
    n_clin = sum(1 for k, v in feature_metadata.items() if v['type'] == 'clinical')
    n_gen  = sum(1 for k, v in feature_metadata.items() if v['type'] == 'genomic')
    print(f"\n  Combined matrix: {len(combined_df)} individuals "
          f"× {combined_df.shape[1]} features "
          f"({n_prs} PRS anchors + {n_clin} clinical + {n_gen} genomic)")

    return combined_df, feature_metadata, per_cohort_iids


def _variance_sparsity_filter(
    df,
    min_variance=0.005,
    max_zero_fraction=0.80,
    max_features=500,
):
    """
    Remove low-information columns before NMF to prevent cluster collapse.

    Drops columns that are:
    - Near-constant after scaling  (variance < min_variance)
    - Predominantly zero           (> max_zero_fraction of values ≤ 1e-9)
    - Beyond top max_features      (by variance, applied last)

    Returns
    -------
    pd.DataFrame  filtered DataFrame
    list          names of retained columns
    """
    n_before = df.shape[1]

    # Near-constant
    var = df.var()
    keep = var[var >= min_variance].index
    df = df[keep]

    # Predominantly zero
    zero_frac = (df <= 1e-9).mean()
    keep = zero_frac[zero_frac <= max_zero_fraction].index
    df = df[keep]

    # Cap by variance (top max_features)
    if df.shape[1] > max_features:
        top_cols = df.var().nlargest(max_features).index
        df = df[top_cols]

    n_after = df.shape[1]
    if n_after < n_before:
        print(f"    Feature filter: {n_before} → {n_after} columns "
              f"(variance≥{min_variance}, sparsity≤{max_zero_fraction*100:.0f}%, "
              f"cap={max_features})")

    return df, df.columns.tolist()


def compute_cohort_cluster_affinity(assignments_df, membership_flags_df, cohorts):
    """
    For each cohort, compute the proportion of its high-risk individuals
    that fall in each cluster.

    Returns
    -------
    pd.DataFrame  (cohorts as index, cluster labels as columns)
        Values in [0, 1] — proportion of cohort's high-risk individuals per cluster.
    """
    cluster_labels = sorted(assignments_df['cluster_label'].unique())
    rows = []

    for cohort in cohorts:
        hr_col = f'is_high_risk_{cohort}'
        if hr_col not in membership_flags_df.columns:
            continue

        hr_iids = membership_flags_df.index[membership_flags_df[hr_col]]
        shared  = assignments_df.index.intersection(hr_iids)

        if len(shared) == 0:
            rows.append({'cohort': cohort, **{c: 0.0 for c in cluster_labels}})
            continue

        counts = assignments_df.loc[shared, 'cluster_label'].value_counts()
        total  = len(shared)
        rows.append({
            'cohort': cohort,
            **{c: counts.get(c, 0) / total for c in cluster_labels},
        })

    affinity_df = pd.DataFrame(rows).set_index('cohort')
    affinity_df = affinity_df[[c for c in cluster_labels if c in affinity_df.columns]]
    return affinity_df


def compute_weighted_feature_list(H_df, feature_metadata):
    """
    For each cluster × feature, compute a weighted importance score:
        weighted_score = H_loading × cohort_weight

    where cohort_weight is |effect_size_r| (clinical) or matching_z_score (genomic).

    Returns
    -------
    pd.DataFrame  tall format with columns:
        cluster, feature, cohort, type, raw_loading, cohort_weight, weighted_score
    Sorted by cluster then weighted_score descending.
    """
    rows = []
    for cluster_label, loadings in H_df.iterrows():
        for feature, loading in loadings.items():
            meta   = feature_metadata.get(feature, {})
            weight = meta.get('weight', 1.0)
            rows.append({
                'cluster':        cluster_label,
                'feature':        feature,
                'cohort':         meta.get('cohort', ''),
                'type':           meta.get('type', ''),
                'raw_loading':    round(float(loading), 6),
                'cohort_weight':  round(float(weight), 6),
                'weighted_score': round(float(loading * weight), 6),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(['cluster', 'weighted_score'], ascending=[True, False])


def annotate_cluster_dominant_cohort(affinity_df):
    """
    For each cluster, identify the cohort with the highest proportion of its
    high-risk individuals in that cluster.

    Returns
    -------
    pd.DataFrame  columns: cluster, dominant_cohort, proportion
    """
    rows = []
    for cluster in affinity_df.columns:
        col = affinity_df[cluster]
        if col.max() == 0:
            rows.append({'cluster': cluster, 'dominant_cohort': None,
                         'proportion': 0.0})
        else:
            dom = col.idxmax()
            rows.append({'cluster': cluster, 'dominant_cohort': dom,
                         'proportion': round(float(col[dom]), 4)})
    return pd.DataFrame(rows)


def plot_cohort_cluster_affinity(affinity_df, output_path):
    """
    Heatmap of cohort × cluster affinity proportions.
    Rows = cohorts, columns = clusters, values = fraction of cohort's
    high-risk individuals in each cluster.
    """
    if affinity_df.empty:
        return

    fig_h = max(3, len(affinity_df) * 0.6)
    fig_w = max(4, len(affinity_df.columns) * 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        affinity_df,
        cmap='YlOrRd',
        vmin=0, vmax=1,
        annot=True, fmt='.2f', annot_kws={'size': 9},
        linewidths=0.5,
        ax=ax,
        cbar_kws={'label': 'Proportion of high-risk individuals'},
    )
    ax.set_title('Cohort–Cluster Affinity\n'
                 '(fraction of each cohort\'s high-risk individuals per cluster)',
                 fontsize=10)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Cohort')
    plt.tight_layout()
    _save_fig(fig, output_path)


# ============================================================================
# PER-COHORT bNMF
# ============================================================================

def run_cohort_bnmf(
    cohort,
    feature_matrix,
    prs_bins,
    k_min=2,
    k_max=6,
    n_runs=30,
    k_select_method='elbow',
    l1_ratio=0.1,
    alpha_W=0.1,
    alpha_H=0.1,
    output_path=None,
    fig_path=None,
    min_variance=0.005,
    max_zero_fraction=0.80,
    max_features=500,
    feature_weights=None,
    scale_method='quantile',
):
    """
    Run consensus bNMF for one cohort and write all outputs.

    Parameters
    ----------
    alpha_H          : float — regularisation strength on H (feature loadings).
                       Must be > 0 when features include sparse/binary-like columns
                       (e.g. HLA allele contributions).  alpha_H=0 allows H to grow
                       unboundedly while W collapses to near-zero, causing everyone
                       to be assigned to cluster_1.  Default 0.1 matches alpha_W.
    min_variance     : float — drop columns with variance below this after scaling
    max_zero_fraction: float — drop columns where more than this fraction of values
                       are zero after scaling (default 0.80 — removes variants
                       present in fewer than 20% of individuals)
    max_features     : int   — cap total features by variance
    feature_weights  : dict | None — {column_name: weight} applied to each feature
                       AFTER impute/scale and variance filtering so weights genuinely
                       influence the NMF rather than being washed out by MinMaxScaler.
                       High-weight features get a proportionally larger range in the
                       NMF input matrix.  None (default) = equal weighting.

    Returns
    -------
    dict  {'assignments', 'H_df', 'coph_df', 'profile_df',
           'optimal_k', 'n'}
    or None if the cohort is skipped / failed.
    """
    if feature_matrix.empty or len(feature_matrix) < MIN_INDIVIDUALS:
        print(f"  [SKIP] {cohort}: only {len(feature_matrix)} individuals "
              f"(minimum {MIN_INDIVIDUALS})")
        return None

    # ── Preprocess ─────────────────────────────────────────────────────────
    X_df, kept_cols = bnmf_core.impute_and_scale(feature_matrix,
                                                   scale_method=scale_method)

    # Remove near-constant and highly sparse columns before NMF
    X_df, kept_cols = _variance_sparsity_filter(
        X_df,
        min_variance=min_variance,
        max_zero_fraction=max_zero_fraction,
        max_features=max_features,
    )
    # ── Post-scale feature weighting ────────────────────────────────────────
    # Applied AFTER MinMaxScaler so weights genuinely stretch the range of
    # high-importance features rather than being normalised away.
    # A feature with weight=5 will have values in [0, 5] instead of [0, 1],
    # giving it 5× more influence in the NMF reconstruction objective.
    if feature_weights:
        applied = 0
        for col in X_df.columns:
            w = feature_weights.get(col, 1.0)
            if w != 1.0 and w > 0:
                X_df[col] = X_df[col] * w
                applied += 1
        if applied:
            print(f"    Weights applied: {applied} features scaled post-normalisation")

    X = X_df.values

    if X.shape[1] < 2:
        print(f"  [SKIP] {cohort}: fewer than 2 features after preprocessing")
        return None

    print(f"\n  Consensus bNMF: {X.shape[0]:,} individuals × {X.shape[1]} features")

    # ── PCA structure diagnostic ────────────────────────────────────────────
    # Run a quick truncated PCA before NMF to assess intrinsic dimensionality.
    # If PC1 explains >70-80% of variance the feature matrix is near rank-1
    # and NMF will always collapse regardless of parameters.  The cumulative
    # variance curve shows how many components are needed to capture the
    # signal and gives an upper bound on meaningful k values for NMF.
    try:
        from sklearn.decomposition import TruncatedSVD
        _n_pca = min(20, X.shape[1] - 1, X.shape[0] - 1)
        # Centre X before SVD — QuantileTransformer gives each feature mean ~0.5,
        # so uncentred SVD would have PC1 capture the mean direction (~0.5 for
        # everyone) rather than actual variance between individuals.
        _X_c = X - X.mean(axis=0)
        _svd = TruncatedSVD(n_components=_n_pca, random_state=42)
        _svd.fit(_X_c)
        _ev = _svd.explained_variance_ratio_
        _cumev = np.cumsum(_ev)
        print(f"  PCA structure (top {_n_pca} components):")
        print(f"    PC1={_ev[0]*100:.1f}%  PC2={_ev[1]*100:.1f}%  "
              f"PC3={_ev[2]*100:.1f}%  PC4={_ev[3]*100:.1f}%  "
              f"PC5={_ev[4]*100:.1f}%")
        for thresh in (0.50, 0.70, 0.90):
            idx = int(np.searchsorted(_cumev, thresh))
            if idx >= len(_cumev):
                print(f"    Components to explain {thresh*100:.0f}% variance: "
                      f">={len(_cumev)} (not reached — re-run with larger --n_pca_components; "
                      f"current cumulative variance = {_cumev[-1]*100:.1f}%)")
            else:
                print(f"    Components to explain {thresh*100:.0f}% variance: {idx + 1}")
        if _ev[0] > 0.70:
            print(f"    *** WARNING: PC1 explains {_ev[0]*100:.1f}% of variance — "
                  f"feature matrix is near rank-1. NMF cluster separation will be "
                  f"poor regardless of k or regularisation. Consider using only "
                  f"clinical features or checking feature diversity.")
        elif _ev[0] > 0.40:
            print(f"    NOTE: PC1 explains {_ev[0]*100:.1f}% — moderate dominance. "
                  f"NMF may find {min(k_max, n_needed)}-{min(k_max, n_needed+2)} "
                  f"meaningful clusters.")
        else:
            print(f"    Good multi-dimensional structure detected — "
                  f"NMF should find meaningful clusters.")
        del _X_c, _svd, _ev, _cumev
    except Exception as _pca_exc:
        print(f"  [WARN] PCA diagnostic failed: {_pca_exc}")

    # Cap k_max to avoid degenerate factorisations
    k_max_eff = min(k_max, X.shape[0] // 5, X.shape[1] - 1)
    k_max_eff = max(k_max_eff, k_min)
    if k_max_eff < k_max:
        print(f"  k_max capped at {k_max_eff} (was {k_max}) "
              f"based on dataset size")

    # ── Consensus NMF ──────────────────────────────────────────────────────
    k_results, coph_df = bnmf_core.run_consensus_bnmf(
        X,
        k_min=k_min, k_max=k_max_eff,
        n_runs=n_runs,
        l1_ratio=l1_ratio,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
    )

    if not k_results or coph_df.empty:
        print(f"  [FAIL] bNMF produced no results for '{cohort}'")
        return None

    # ── Select optimal k ───────────────────────────────────────────────────
    optimal_k = bnmf_core.select_optimal_k(coph_df, method=k_select_method)

    # ── Cluster assignments + feature loadings ─────────────────────────────
    assignments_df, H_df = bnmf_core.get_consensus_cluster_assignments(
        k_results, optimal_k,
        feature_names=kept_cols,
        index=X_df.index,
    )

    # Attach phenotype labels
    if 'PHENOTYPE' in prs_bins.columns:
        assignments_df = assignments_df.join(
            prs_bins[['PHENOTYPE']], how='left'
        )

    # ── Cluster profiles ───────────────────────────────────────────────────
    phenotype_series = (prs_bins['PHENOTYPE']
                        if 'PHENOTYPE' in prs_bins.columns else None)
    profile_df = bnmf_core.characterize_clusters(
        assignments_df, X_df,
        phenotype_series=phenotype_series,
    )

    # ── Save outputs ───────────────────────────────────────────────────────
    if output_path:
        os.makedirs(output_path, exist_ok=True)

        assignments_df.to_csv(
            os.path.join(output_path, 'cluster_assignments.csv'))
        H_df.to_csv(
            os.path.join(output_path, 'feature_loadings.csv'))
        coph_df.to_csv(
            os.path.join(output_path, 'cophenetic_curve.csv'), index=False)

        # Basis matrix: normalised membership columns from assignments
        membership_cols = [c for c in assignments_df.columns
                           if c.startswith('membership_')]
        assignments_df[membership_cols].to_csv(
            os.path.join(output_path, 'basis_matrix.csv'))

        if not profile_df.empty:
            profile_df.to_csv(
                os.path.join(output_path, 'cluster_profile.csv'), index=True)  # index='cluster'

        print(f"    Outputs → {output_path}")

    # ── Figures ────────────────────────────────────────────────────────────
    if fig_path:
        os.makedirs(fig_path, exist_ok=True)

        plot_cophenetic_curve(
            coph_df,
            os.path.join(fig_path, 'cophenetic_curve.png'),
            cohort,
        )
        # Recompute C from W_list — it is not stored in k_results to save memory
        C_optimal = bnmf_core._soft_consensus_matrix(k_results[optimal_k]['W_list'])
        plot_consensus_map(
            C_optimal, optimal_k,
            os.path.join(fig_path, f'consensus_map_k{optimal_k}.png'),
            cohort,
        )
        del C_optimal
        if not profile_df.empty:
            plot_cluster_profile(
                profile_df,
                os.path.join(fig_path, 'cluster_profile.png'),
                cohort,
            )
        plot_feature_loadings(
            H_df,
            os.path.join(fig_path, 'feature_loadings.png'),
            cohort,
        )

    return {
        'assignments': assignments_df,
        'H_df':        H_df,
        'coph_df':     coph_df,
        'profile_df':  profile_df,
        'optimal_k':   optimal_k,
        'n':           len(assignments_df),
    }


# ============================================================================
# REPLOT FROM SAVED CSVs
# ============================================================================

def replot_figures_from_saved_csvs(pheno_data, top_n=20):
    """
    Regenerate cluster_profile.png and feature_loadings.png for every
    result directory that already has the corresponding CSV files.

    Scans both per-cohort directories:
        {pheno_data}/scores/cohortBnmf/{cohort}/{population}/
    and combined directories:
        {pheno_data}/scores/combinedCohortBnmf_{population}/

    No NMF is re-run — only plotting functions are called.
    """
    scores_root = os.path.join(pheno_data, 'scores')
    figs_root   = os.path.join(pheno_data, 'figures')
    n_replotted = 0

    def _replot_dir(result_dir, fig_dir, cohort_label):
        nonlocal n_replotted
        profile_csv  = os.path.join(result_dir, 'cluster_profile.csv')
        loadings_csv = os.path.join(result_dir, 'feature_loadings.csv')
        os.makedirs(fig_dir, exist_ok=True)

        if os.path.exists(profile_csv):
            try:
                profile_df = pd.read_csv(profile_csv, index_col='cluster')
                plot_cluster_profile(
                    profile_df,
                    os.path.join(fig_dir, 'cluster_profile.png'),
                    cohort_label,
                    top_n=top_n,
                )
                n_replotted += 1
            except Exception as exc:
                print(f"  [WARN] cluster_profile replot failed for {cohort_label}: {exc}")

        if os.path.exists(loadings_csv):
            try:
                H_df = pd.read_csv(loadings_csv, index_col=0)
                plot_feature_loadings(
                    H_df,
                    os.path.join(fig_dir, 'feature_loadings.png'),
                    cohort_label,
                    top_n=top_n,
                )
                n_replotted += 1
            except Exception as exc:
                print(f"  [WARN] feature_loadings replot failed for {cohort_label}: {exc}")

    # ── Per-cohort directories ───────────────────────────────────────────────
    cohort_bnmf_root = os.path.join(scores_root, 'cohortBnmf')
    if os.path.isdir(cohort_bnmf_root):
        for cohort in sorted(os.listdir(cohort_bnmf_root)):
            cohort_dir = os.path.join(cohort_bnmf_root, cohort)
            if not os.path.isdir(cohort_dir):
                continue
            for population in sorted(os.listdir(cohort_dir)):
                pop_dir = os.path.join(cohort_dir, population)
                if not os.path.isdir(pop_dir):
                    continue
                fig_dir = os.path.join(figs_root, 'cohortBnmf', cohort, population)
                label   = f'{cohort}/{population}'
                print(f"  Replotting {label} ...")
                _replot_dir(pop_dir, fig_dir, label)

    # ── Combined directories ─────────────────────────────────────────────────
    for entry in sorted(os.listdir(scores_root)):
        if not entry.startswith('combinedCohortBnmf'):
            continue
        combined_dir = os.path.join(scores_root, entry)
        if not os.path.isdir(combined_dir):
            continue
        fig_dir = os.path.join(figs_root, entry)
        print(f"  Replotting {entry} ...")
        _replot_dir(combined_dir, fig_dir, entry)

    print(f"\nReplot complete — {n_replotted} figure(s) regenerated.")


# ============================================================================
# MAIN ORCHESTRATOR
# ============================================================================

def run_all_cohorts(
    pheno_data,
    raw_features_file,
    filter_strategy='tiered',
    use_set='validation',
    k_min=2,
    k_max=6,
    n_runs=30,
    k_select_method='elbow',
    l1_ratio=0.1,
    alpha_W=0.1,
    alpha_H=0.1,
    min_effect_size=0.0,
    specificity_tier_max=3,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    cohorts_to_run=None,
    include_holdout_prs=True,
    population='both',
    min_variance=0.005,
    max_zero_fraction=0.80,
    max_features=500,
    include_raw_clinical=True,
    scale_method='quantile',
):
    """
    Run cohort bNMF for all (or specified) cohorts.

    Parameters
    ----------
    population : str  'both' | 'high_risk' | 'low_control' | 'separate'
        Which population to include in the feature matrix:
          'both'        — union of high-risk cases and low-control controls
                          (original behaviour, single output per cohort)
          'high_risk'   — only high-risk cases (bin > threshold, PHENOTYPE=case)
          'low_control' — only low-risk controls (bin < threshold, PHENOTYPE=ctrl)
          'separate'    — run two independent bNMF analyses per cohort; saves
                          results to {cohort}/high_risk/ and {cohort}/low_control/
                          so each population gets its own clusters and enrichment.

    Returns
    -------
    pd.DataFrame  overview table (one row per cohort × population)
    """
    scores_path  = os.path.join(pheno_data, 'scores')
    base_output  = os.path.join(scores_path, 'cohortBnmf')
    base_fig     = os.path.join(pheno_data,  'figures', 'cohortBnmf')
    os.makedirs(base_output, exist_ok=True)
    os.makedirs(base_fig,    exist_ok=True)

    print("=" * 60)
    print("COHORT bNMF SUBTYPE ANALYSIS")
    print("=" * 60)
    print(f"  pheno_data            : {pheno_data}")
    print(f"  raw_features_file     : {raw_features_file}")
    print(f"  filter_strategy       : {filter_strategy}")
    print(f"  use_set               : {use_set}")
    print(f"  k_min / k_max         : {k_min} / {k_max}")
    print(f"  n_runs                : {n_runs}")
    print(f"  k_select_method       : {k_select_method}")
    print(f"  high_risk_threshold   : bin > {high_risk_threshold}")
    print(f"  low_control_threshold : bin < {low_control_threshold}")
    print(f"  min_effect_size       : {min_effect_size}")
    print(f"  specificity_tier_max  : {specificity_tier_max}")
    print(f"  population            : {population}")

    # ── Load shared data ───────────────────────────────────────────────────
    print("\n[1] Loading cohort analysis results...")
    genomic_df  = load_cohort_genomic_features(scores_path, filter_strategy)
    clinical_df = load_cohort_clinical_features(scores_path, use_set)

    print("\n[2] Loading individual-level data...")
    prs_bins     = load_prs_bins(scores_path, use_set)
    raw_clinical = load_raw_clinical_data(raw_features_file)
    prs_dict     = load_all_prs_feature_files(
        scores_path, include_holdout=include_holdout_prs, use_set=use_set)

    # ── Determine cohorts ──────────────────────────────────────────────────
    bin_cols   = [c for c in prs_bins.columns
                  if c.startswith('bin_') and 'combined' not in c]
    all_cohorts = [c.replace('bin_', '') for c in bin_cols]

    if cohorts_to_run:
        cohorts = [c for c in all_cohorts if c in cohorts_to_run]
        unknown = [c for c in cohorts_to_run if c not in all_cohorts]
        if unknown:
            print(f"  [WARN] Requested cohorts not found in PRS bins: {unknown}")
    else:
        cohorts = all_cohorts

    print(f"\n[3] Cohorts to analyse: {cohorts}")

    # Filter prs_dict to models matching valid cohorts only.
    # Removes covariates, all_*, all+env_combined, etc.
    prs_dict = _filter_prs_dict_to_cohorts(prs_dict, cohorts)

    # ── Per-cohort loop ────────────────────────────────────────────────────
    # When population='separate', each cohort is run twice (high_risk then
    # low_control) with independent bNMF analyses saved to separate subdirs.
    overview_rows = []
    pop_runs = (['high_risk', 'low_control']
                if population == 'separate' else [population])

    for cohort in cohorts:
        for pop in pop_runs:
            pop_label = f" [{pop}]" if population == 'separate' else ''
            print(f"\n{'='*60}")
            print(f"  COHORT: {cohort}{pop_label}")
            print(f"{'='*60}")

            # Individual selection
            selected_iids, sel_stats = select_cohort_individuals(
                prs_bins, cohort,
                high_risk_threshold=high_risk_threshold,
                low_control_threshold=low_control_threshold,
                population=pop,
            )

            if len(selected_iids) < MIN_INDIVIDUALS:
                print(f"  [SKIP] Too few individuals: "
                      f"{len(selected_iids)} < {MIN_INDIVIDUALS}")
                overview_rows.append({
                    'cohort':     cohort,
                    'population': pop,
                    'status':     'skipped_too_few_individuals',
                    **sel_stats,
                })
                continue

            # Feature matrix
            print(f"  Building feature matrix...")
            feature_matrix, feat_counts = build_cohort_feature_matrix(
                cohort=cohort,
                selected_iids=selected_iids,
                genomic_features_df=genomic_df,
                clinical_features_df=clinical_df,
                prs_dict=prs_dict,
                raw_clinical_df=raw_clinical,
                min_effect_size=min_effect_size,
                specificity_tier_max=specificity_tier_max,
                include_raw_clinical=include_raw_clinical,
            )

            if feature_matrix.empty or feature_matrix.shape[1] < 2:
                print(f"  [SKIP] Insufficient features for '{cohort}'")
                overview_rows.append({
                    'cohort':     cohort,
                    'population': pop,
                    'status':     'skipped_no_features',
                    **sel_stats, **feat_counts,
                })
                continue

            # Output paths — subdirectory per population when running separately
            if population == 'separate':
                cohort_output = os.path.join(base_output, cohort, pop)
                cohort_figs   = os.path.join(base_fig,    cohort, pop)
            else:
                cohort_output = os.path.join(base_output, cohort)
                cohort_figs   = os.path.join(base_fig,    cohort)

            result = run_cohort_bnmf(
                cohort=cohort,
                feature_matrix=feature_matrix,
                prs_bins=prs_bins,
                k_min=k_min,
                k_max=k_max,
                n_runs=n_runs,
                k_select_method=k_select_method,
                l1_ratio=l1_ratio,
                alpha_W=alpha_W,
                alpha_H=alpha_H,
                output_path=cohort_output,
                fig_path=cohort_figs,
                min_variance=min_variance,
                max_zero_fraction=max_zero_fraction,
                max_features=max_features,
                scale_method=scale_method,
            )

            if result:
                print(f"  Completed: k={result['optimal_k']}, n={result['n']:,}")
                overview_rows.append({
                    'cohort':         cohort,
                    'population':     pop,
                    'status':         'completed',
                    'optimal_k':      result['optimal_k'],
                    'n_in_analysis':  result['n'],
                    **sel_stats,
                    **feat_counts,
                })
            else:
                overview_rows.append({
                    'cohort':     cohort,
                    'population': pop,
                    'status':     'failed',
                    **sel_stats,
                    **feat_counts,
                })

    # ── Overview table ─────────────────────────────────────────────────────
    overview_df = pd.DataFrame(overview_rows)
    overview_path = os.path.join(base_output, 'cohort_bnmf_overview.csv')
    overview_df.to_csv(overview_path, index=False)

    print("\n" + "=" * 60)
    print("COHORT bNMF COMPLETE")
    print("=" * 60)

    display_cols = ['cohort', 'status', 'n_total']
    if 'optimal_k' in overview_df.columns:
        display_cols.append('optimal_k')
    if 'n_in_analysis' in overview_df.columns:
        display_cols.append('n_in_analysis')
    print(overview_df[display_cols].to_string(index=False))
    print(f"\nOverview → {overview_path}")
    print(f"Outputs  → {base_output}/{{cohort}}/")
    print(f"Figures  → {base_fig}/{{cohort}}/")

    return overview_df


# ============================================================================
# PCA STRUCTURE ANALYSIS
# ============================================================================

def run_pca_analysis(
    pheno_data,
    raw_features_file,
    filter_strategy='tiered',
    use_set='validation',
    min_effect_size=0.0,
    specificity_tier_max=3,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    cohorts_to_run=None,
    include_holdout_prs=True,
    min_variance=0.005,
    max_zero_fraction=0.80,
    max_features=500,
    include_prs_with_genomic=False,
    population='high_risk',
    include_raw_clinical=True,
    scale_method='quantile',
    n_pca_components=30,
    clinical_only=False,
):
    """
    Build the same combined feature matrix used by run_combined_bnmf, scale
    it identically, then run PCA to assess intrinsic dimensionality before
    committing to NMF.

    Outputs (scores/combinedPCA_{population}/)
    -------------------------------------------
    pca_scores.csv        — individual × PC coordinates (PC1..PCn)
    pca_loadings.csv      — feature × PC component weights
    pca_variance.csv      — explained variance ratio per PC + cumulative

    Figures (figures/combinedPCA_{population}/)
    --------------------------------------------
    pca_scree.png         — per-PC explained variance + cumulative curve
    pca_scatter_cohort.png — PC1 vs PC2 coloured by cohort membership
    pca_scatter_prs.png   — PC1 vs PC2 coloured by combined PRS bin
    pca_loadings.png      — top feature loadings for PC1–PC4
    """
    from sklearn.decomposition import TruncatedSVD

    scores_path = os.path.join(pheno_data, 'scores')
    pop_suffix  = f'_{population}' if population != 'both' else ''
    feat_suffix = '_clinical_only' if clinical_only else ''
    base_output = os.path.join(scores_path,           f'combinedPCA{pop_suffix}{feat_suffix}')
    base_fig    = os.path.join(pheno_data, 'figures', f'combinedPCA{pop_suffix}{feat_suffix}')
    os.makedirs(base_output, exist_ok=True)
    os.makedirs(base_fig,    exist_ok=True)

    print("\n" + "=" * 60)
    print("COMBINED PCA STRUCTURE ANALYSIS")
    print("=" * 60)
    print(f"  population    : {population}")
    print(f"  use_set       : {use_set}")
    print(f"  features      : {'clinical only' if clinical_only else 'genomic + clinical'}")

    # ── Load shared data (identical to run_combined_bnmf) ─────────────────
    print("\n[1] Loading feature definitions...")
    genomic_df  = load_cohort_genomic_features(scores_path, filter_strategy)
    clinical_df = load_cohort_clinical_features(scores_path, use_set)

    print("\n[2] Loading individual-level data...")
    prs_bins     = load_prs_bins(scores_path, use_set)
    raw_clinical = load_raw_clinical_data(raw_features_file)
    prs_dict     = load_all_prs_feature_files(
        scores_path, include_holdout=include_holdout_prs, use_set=use_set)

    bin_cols    = [c for c in prs_bins.columns
                   if c.startswith('bin_') and 'combined' not in c]
    all_cohorts = [c.replace('bin_', '') for c in bin_cols]
    cohorts     = ([c for c in all_cohorts if c in cohorts_to_run]
                   if cohorts_to_run else all_cohorts)
    prs_dict    = _filter_prs_dict_to_cohorts(prs_dict, cohorts)
    print(f"\n[3] Cohorts: {cohorts}")

    # ── Build combined feature matrix ─────────────────────────────────────
    print("\n[4] Building combined feature matrix...")
    combined_matrix, _, _ = build_deduplicated_feature_matrix(
        cohorts=cohorts,
        prs_bins=prs_bins,
        genomic_features_df=genomic_df,
        clinical_features_df=clinical_df,
        prs_dict=prs_dict,
        raw_clinical_df=raw_clinical,
        high_risk_threshold=high_risk_threshold,
        low_control_threshold=low_control_threshold,
        min_effect_size=min_effect_size,
        specificity_tier_max=specificity_tier_max,
        include_prs_with_genomic=include_prs_with_genomic,
        population=population,
        include_raw_clinical=include_raw_clinical,
    )

    if combined_matrix.empty:
        print("  [ERROR] Combined feature matrix is empty. Aborting.")
        return

    # ── Clinical-only filter ───────────────────────────────────────────────
    if clinical_only:
        clin_cols = [c for c in combined_matrix.columns if c.startswith('clin_')]
        if not clin_cols:
            print("  [ERROR] No clinical (clin_*) columns found in combined matrix. "
                  "Ensure --no_raw_clinical is NOT set when using --clinical_only.")
            return
        combined_matrix = combined_matrix[clin_cols]
        print(f"  Clinical-only: retained {len(clin_cols)} clin_* columns, "
              f"dropped all gen_* and PRS anchor columns.")

    print(f"\n[5] Scaling ({scale_method}) and filtering...")
    X_df, _retained = bnmf_core.impute_and_scale(
        combined_matrix, scale_method=scale_method)
    X_df, _cols = _variance_sparsity_filter(
        X_df, min_variance=min_variance,
        max_zero_fraction=max_zero_fraction, max_features=max_features)
    print(f"  Matrix for PCA: {X_df.shape[0]:,} individuals × {X_df.shape[1]} features")

    # ── Cohort membership flags for colouring ─────────────────────────────
    membership_flags = build_cohort_membership_flags(
        prs_bins, cohorts, high_risk_threshold, low_control_threshold)

    # Assign a primary cohort label per individual (first high-risk cohort found)
    def _primary_cohort(iid):
        if iid not in membership_flags.index:
            return 'none'
        row = membership_flags.loc[iid]
        hr_cols = [c for c in row.index if c.startswith('is_high_risk_')]
        for col in hr_cols:
            if row[col]:
                return col.replace('is_high_risk_', '')
        lc_cols = [c for c in row.index if c.startswith('is_low_control_')]
        for col in lc_cols:
            if row[col]:
                return col.replace('is_low_control_', '') + '_low'
        return 'other'

    cohort_labels = pd.Series(
        [_primary_cohort(iid) for iid in X_df.index],
        index=X_df.index, name='primary_cohort',
    )
    # Count how many cohorts each individual is high-risk in
    hr_flag_cols = [c for c in membership_flags.columns
                    if c.startswith('is_high_risk_')]
    n_cohorts_hr = membership_flags.reindex(X_df.index)[hr_flag_cols].sum(axis=1).fillna(0)

    # ── Run PCA (centred TruncatedSVD) ────────────────────────────────────
    print("\n[6] Running PCA...")
    X = X_df.values.astype(np.float32)
    X_c = X - X.mean(axis=0)
    n_components = min(n_pca_components, X.shape[1] - 1, X.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=42)
    scores = svd.fit_transform(X_c)       # n_individuals × n_components
    loadings = svd.components_            # n_components × n_features
    ev  = svd.explained_variance_ratio_
    cev = np.cumsum(ev)

    print(f"  Explained variance:")
    print(f"    PC1={ev[0]*100:.1f}%  PC2={ev[1]*100:.1f}%  "
          f"PC3={ev[2]*100:.1f}%  PC4={ev[3]*100:.1f}%  PC5={ev[4]*100:.1f}%")
    for thresh in (0.50, 0.70, 0.90):
        idx = int(np.searchsorted(cev, thresh))
        if idx >= len(cev):
            print(f"    PCs to explain {thresh*100:.0f}%: >={len(cev)} "
                  f"(not reached — current cumulative = {cev[-1]*100:.1f}%; "
                  f"re-run with larger --n_pca_components)")
        else:
            print(f"    PCs to explain {thresh*100:.0f}%: {idx + 1}")
    if ev[0] > 0.50:
        print(f"  *** PC1 explains {ev[0]*100:.1f}% — strong single dominant direction. "
              f"NMF cluster separation will be limited.")
    elif ev[0] > 0.25:
        print(f"  NOTE: PC1 explains {ev[0]*100:.1f}% — moderate. Check scatter plots "
              f"for visible sub-clusters.")
    else:
        print(f"  Good multi-dimensional structure — NMF should find meaningful clusters.")

    # ── Save CSVs ──────────────────────────────────────────────────────────
    pc_cols = [f'PC{i+1}' for i in range(n_components)]
    scores_df = pd.DataFrame(scores, index=X_df.index, columns=pc_cols)
    scores_df['primary_cohort'] = cohort_labels
    scores_df['n_cohorts_high_risk'] = n_cohorts_hr
    scores_df.to_csv(os.path.join(base_output, 'pca_scores.csv'))

    loadings_df = pd.DataFrame(loadings, index=pc_cols, columns=X_df.columns)
    loadings_df.to_csv(os.path.join(base_output, 'pca_loadings.csv'))

    variance_df = pd.DataFrame({
        'PC': pc_cols,
        'explained_variance_ratio': ev,
        'cumulative_variance': cev,
    })
    variance_df.to_csv(os.path.join(base_output, 'pca_variance.csv'), index=False)
    print(f"  CSVs → {base_output}")

    # ── Figure 1: Scree + cumulative variance ─────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    pc_nums = list(range(1, n_components + 1))

    axes[0].bar(pc_nums, ev * 100, color='#1f77b4', alpha=0.8)
    axes[0].set_xlabel('Principal Component')
    axes[0].set_ylabel('Explained variance (%)')
    axes[0].set_title('Per-PC Explained Variance')
    axes[0].axhline(ev[0] * 100, color='red', linestyle='--', linewidth=0.8, alpha=0.5)

    axes[1].plot(pc_nums, cev * 100, marker='o', color='#1f77b4',
                 linewidth=2, markersize=4)
    for thresh in (50, 70, 90):
        idx = int(np.searchsorted(cev, thresh / 100))
        axes[1].axhline(thresh, color='grey', linestyle=':', linewidth=0.8)
        if idx < len(cev):
            axes[1].annotate(f'{thresh}% @ PC{idx + 1}',
                             xy=(idx + 1, thresh), xytext=(4, 2),
                             textcoords='offset points', fontsize=7, color='grey')
        else:
            # Threshold not reached — annotate at right edge with ">N" label
            axes[1].annotate(f'{thresh}% not reached\n(>{len(cev)} PCs)',
                             xy=(len(cev), cev[-1] * 100), xytext=(-60, 4),
                             textcoords='offset points', fontsize=6, color='red')
    axes[1].set_xlabel('Number of Components')
    axes[1].set_ylabel('Cumulative variance (%)')
    axes[1].set_title('Cumulative Explained Variance')
    axes[1].set_ylim(0, 105)

    plt.suptitle(f'PCA Structure — combined {population} ({X_df.shape[0]:,} individuals × '
                 f'{X_df.shape[1]} features)', fontsize=10)
    plt.tight_layout()
    _save_fig(fig, os.path.join(base_fig, 'pca_scree.png'))

    # ── Figure 2: PC1 vs PC2 coloured by cohort ───────────────────────────
    unique_cohorts = sorted(cohort_labels.unique())
    palette = plt.cm.get_cmap('tab10', len(unique_cohorts))
    colour_map = {c: palette(i) for i, c in enumerate(unique_cohorts)}

    fig, ax = plt.subplots(figsize=(9, 7))
    for cohort_name in unique_cohorts:
        mask = cohort_labels == cohort_name
        ax.scatter(
            scores[mask, 0], scores[mask, 1],
            c=[colour_map[cohort_name]], label=cohort_name,
            alpha=0.4, s=8, linewidths=0,
        )
    ax.set_xlabel(f'PC1  ({ev[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2  ({ev[1]*100:.1f}%)')
    ax.set_title(f'PC1 vs PC2 — coloured by primary cohort\n'
                 f'({population}, {X_df.shape[0]:,} individuals)')
    ax.legend(fontsize=7, markerscale=3, loc='best')
    plt.tight_layout()
    _save_fig(fig, os.path.join(base_fig, 'pca_scatter_cohort.png'))

    # ── Figure 3: PC1 vs PC2 coloured by number of cohorts (high-risk) ────
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        scores[:, 0], scores[:, 1],
        c=n_cohorts_hr.values, cmap='viridis',
        alpha=0.4, s=8, linewidths=0,
    )
    plt.colorbar(sc, ax=ax, label='# cohorts individual is high-risk in')
    ax.set_xlabel(f'PC1  ({ev[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2  ({ev[1]*100:.1f}%)')
    ax.set_title(f'PC1 vs PC2 — coloured by multi-cohort high-risk membership\n'
                 f'({population}, {X_df.shape[0]:,} individuals)')
    plt.tight_layout()
    _save_fig(fig, os.path.join(base_fig, 'pca_scatter_ncohorts.png'))

    # ── Figure 4: Top feature loadings for PC1–PC4 ────────────────────────
    n_show   = min(4, n_components)
    top_n_ft = 20
    fig, axes = plt.subplots(1, n_show, figsize=(n_show * 5, 7))
    if n_show == 1:
        axes = [axes]
    for i, ax in enumerate(axes[:n_show]):
        pc_load = pd.Series(loadings[i], index=X_df.columns)
        top_pos = pc_load.nlargest(top_n_ft // 2)
        top_neg = pc_load.nsmallest(top_n_ft // 2)
        top     = pd.concat([top_pos, top_neg]).sort_values()
        colours = ['#d62728' if f.startswith('clin_') else '#1f77b4'
                   for f in top.index]
        labels  = [f.replace('clin_', 'C: ').replace('gen_', 'G: ')
                   for f in top.index]
        ax.barh(range(len(top)), top.values, color=colours, alpha=0.85)
        ax.axvline(0, color='black', linewidth=0.6, linestyle='--')
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(f'PC{i+1}  ({ev[i]*100:.1f}%)', fontsize=9)
        ax.set_xlabel('Loading')
    legend_handles = [
        mpatches.Patch(facecolor='#d62728', label='Clinical'),
        mpatches.Patch(facecolor='#1f77b4', label='Genomic'),
    ]
    axes[-1].legend(handles=legend_handles, fontsize=7, loc='lower right')
    plt.suptitle(f'Top Feature Loadings per PC — combined {population}', fontsize=10)
    plt.tight_layout()
    _save_fig(fig, os.path.join(base_fig, 'pca_loadings.png'))

    print(f"\n  Figures → {base_fig}")
    print("\n" + "=" * 60)
    print("PCA COMPLETE")
    print("=" * 60)


# ============================================================================
# COMBINED ORCHESTRATOR
# ============================================================================

def run_combined_bnmf(
    pheno_data,
    raw_features_file,
    filter_strategy='tiered',
    use_set='validation',
    k_min=2,
    k_max=6,
    n_runs=30,
    k_select_method='elbow',
    l1_ratio=0.1,
    alpha_W=0.1,
    alpha_H=0.1,
    min_effect_size=0.0,
    specificity_tier_max=3,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    cohorts_to_run=None,
    include_holdout_prs=True,
    weight_features=False,
    min_variance=0.005,
    max_zero_fraction=0.80,
    max_features=500,
    include_prs_with_genomic=False,
    population='both',
    include_raw_clinical=True,
    scale_method='quantile',
    clinical_only=False,
):
    """
    Combined cross-cohort bNMF subtype analysis.

    Builds a single de-duplicated feature matrix:
      {cohort}_prs   — scaled_prs per cohort model (non-redundant anchors)
      clin_{feature} — raw clinical value, one column; weight = max |effect_size_r|
      gen_{feature}  — weighted contribution from best model; weight = max matching_z_score

    This avoids the cluster-collapse caused by cohort-prefixed redundant columns.

    population : str  'both' | 'high_risk' | 'low_control' | 'separate'
        Controls which individuals are included across all cohorts.
        'separate' runs two combined analyses — one high_risk, one low_control.

    Post-analysis outputs
    ---------------------
    cohort_cluster_affinity.csv  — cohort × cluster proportions of high-risk individuals
    weighted_feature_list.csv    — per-cluster features ranked by loading × cohort weight
    cluster_dominant_cohort.csv  — which cohort dominates each cluster
    cluster_assignments.csv      — includes is_high_risk_{cohort} / is_low_control_{cohort} flags
    feature_weights.csv          — importance weight per feature used for weighting
    cohort_cluster_affinity.png  — heatmap of affinity matrix
    """
    scores_path = os.path.join(pheno_data, 'scores')
    # Use a population-specific subdirectory so high_risk and low_control
    # combined runs don't overwrite each other's outputs.
    pop_suffix  = f'_{population}' if population != 'both' else ''
    clin_suffix = '_clinical_only' if clinical_only else ''
    base_output = os.path.join(scores_path, f'combinedCohortBnmf{pop_suffix}{clin_suffix}')
    base_fig    = os.path.join(pheno_data,  'figures', f'combinedCohortBnmf{pop_suffix}{clin_suffix}')

    print("\n" + "=" * 60)
    print("COMBINED CROSS-COHORT bNMF")
    print("=" * 60)

    # ── Load shared data ───────────────────────────────────────────────────
    print("\n[1] Loading feature definitions...")
    genomic_df  = load_cohort_genomic_features(scores_path, filter_strategy)
    clinical_df = load_cohort_clinical_features(scores_path, use_set)

    print("\n[2] Loading individual-level data...")
    prs_bins     = load_prs_bins(scores_path, use_set)
    raw_clinical = load_raw_clinical_data(raw_features_file)
    prs_dict     = load_all_prs_feature_files(
        scores_path, include_holdout=include_holdout_prs, use_set=use_set)

    # ── Determine cohorts ──────────────────────────────────────────────────
    bin_cols    = [c for c in prs_bins.columns
                   if c.startswith('bin_') and 'combined' not in c]
    all_cohorts = [c.replace('bin_', '') for c in bin_cols]

    if cohorts_to_run:
        cohorts = [c for c in all_cohorts if c in cohorts_to_run]
        unknown = [c for c in cohorts_to_run if c not in all_cohorts]
        if unknown:
            print(f"  [WARN] Requested cohorts not found in PRS bins: {unknown}")
    else:
        cohorts = all_cohorts

    prs_dict = _filter_prs_dict_to_cohorts(prs_dict, cohorts)
    print(f"\n[3] Cohorts contributing to combined matrix: {cohorts}")

    # ── Handle 'separate' population: recurse for each sub-population ──────
    if population == 'separate':
        print("\n  population='separate' → running combined bNMF twice "
              "(high_risk then low_control)")
        kwargs = dict(
            pheno_data=pheno_data, raw_features_file=raw_features_file,
            filter_strategy=filter_strategy, use_set=use_set,
            k_min=k_min, k_max=k_max, n_runs=n_runs,
            k_select_method=k_select_method, l1_ratio=l1_ratio,
            alpha_W=alpha_W, alpha_H=alpha_H,
            min_effect_size=min_effect_size,
            specificity_tier_max=specificity_tier_max,
            high_risk_threshold=high_risk_threshold,
            low_control_threshold=low_control_threshold,
            cohorts_to_run=cohorts_to_run,
            include_holdout_prs=include_holdout_prs,
            weight_features=weight_features,
            min_variance=min_variance, max_zero_fraction=max_zero_fraction,
            max_features=max_features,
            include_prs_with_genomic=include_prs_with_genomic,
            include_raw_clinical=include_raw_clinical,
            scale_method=scale_method,
        )
        run_combined_bnmf(**kwargs, population='high_risk')
        run_combined_bnmf(**kwargs, population='low_control')
        return

    # ── Build combined feature matrix (de-duplicated) ─────────────────────
    print("\n[4] Building de-duplicated feature matrix...")
    print(f"    Population filter: '{population}'")
    combined_matrix, feature_metadata, per_cohort_iids = build_deduplicated_feature_matrix(
        cohorts=cohorts,
        prs_bins=prs_bins,
        genomic_features_df=genomic_df,
        clinical_features_df=clinical_df,
        prs_dict=prs_dict,
        raw_clinical_df=raw_clinical,
        high_risk_threshold=high_risk_threshold,
        low_control_threshold=low_control_threshold,
        min_effect_size=min_effect_size,
        specificity_tier_max=specificity_tier_max,
        include_prs_with_genomic=include_prs_with_genomic,
        population=population,
        include_raw_clinical=include_raw_clinical,
    )

    if combined_matrix.empty:
        print("  [ERROR] Combined feature matrix is empty. Aborting.")
        return None

    if len(combined_matrix) < MIN_INDIVIDUALS:
        print(f"  [ERROR] Too few individuals ({len(combined_matrix)} < {MIN_INDIVIDUALS}). Aborting.")
        return None

    # ── Build cohort membership flags ─────────────────────────────────────
    membership_flags = build_cohort_membership_flags(
        prs_bins, cohorts, high_risk_threshold, low_control_threshold
    )

    # ── Feature weighting (post-scale, inside run_cohort_bnmf) ────────────
    # Weights are passed to run_cohort_bnmf() and applied AFTER MinMaxScaler
    # so they genuinely stretch high-importance feature ranges rather than
    # being normalised away.  When weight_features=False (default), None is
    # passed and all features are treated equally.
    nmf_matrix = combined_matrix.copy()

    # ── Clinical-only mode: drop all genomic features ──────────────────────
    if clinical_only:
        clin_cols = [c for c in nmf_matrix.columns if c.startswith('clin_')]
        if not clin_cols:
            print("  [ERROR] No clinical (clin_*) columns found. "
                  "Ensure --no_raw_clinical is NOT set when using --clinical_only.")
            return
        dropped_gen = len(nmf_matrix.columns) - len(clin_cols)
        nmf_matrix = nmf_matrix[clin_cols]
        print(f"  Clinical-only: {len(clin_cols)} clin_* features retained, "
              f"{dropped_gen} gen_* features dropped")

    resolved_weights = None
    if weight_features and feature_metadata:
        resolved_weights = {
            col: feature_metadata[col].get('weight', 1.0)
            for col in nmf_matrix.columns
            if col in feature_metadata
        }
        print(f"    Feature importance weighting enabled "
              f"({sum(1 for w in resolved_weights.values() if w != 1.0)} "
              f"features with non-unit weight) — applied post-scaling inside NMF")

    # Save weight metadata regardless (useful for inspection even when not applied)
    if feature_metadata:
        os.makedirs(base_output, exist_ok=True)
        weight_df = pd.DataFrame([
            {'feature': col, **feature_metadata.get(col, {})}
            for col in nmf_matrix.columns
        ])
        weight_df.to_csv(os.path.join(base_output, 'feature_weights.csv'), index=False)

    # ── Run bNMF on combined matrix ────────────────────────────────────────
    print("\n[5] Running bNMF on combined matrix...")
    print(f"    Matrix shape before impute/scale: {nmf_matrix.shape}")
    print(f"    Variance/sparsity filter params: "
          f"min_var={min_variance}, max_zero={max_zero_fraction}, cap={max_features}")
    os.makedirs(base_output, exist_ok=True)
    os.makedirs(base_fig,    exist_ok=True)

    result = run_cohort_bnmf(
        cohort='combined',
        feature_matrix=nmf_matrix,
        prs_bins=prs_bins,
        min_variance=min_variance,
        max_zero_fraction=max_zero_fraction,
        max_features=max_features,
        output_path=base_output,
        fig_path=base_fig,
        k_min=k_min,
        k_max=k_max,
        n_runs=n_runs,
        k_select_method=k_select_method,
        l1_ratio=l1_ratio,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
        feature_weights=resolved_weights,
        scale_method=scale_method,
    )

    if result is None:
        print("  [ERROR] bNMF run failed.")
        return None

    assignments_df = result['assignments']
    H_df           = result['H_df']

    # ── Post-analysis ──────────────────────────────────────────────────────
    print("\n[6] Post-analysis...")

    # (a) Attach cohort membership flags to cluster assignments
    flag_cols = [c for c in membership_flags.columns
                 if c in [f'is_high_risk_{ch}' for ch in cohorts]
                 + [f'is_low_control_{ch}' for ch in cohorts]]
    shared_idx = assignments_df.index.intersection(membership_flags.index)
    if len(shared_idx) > 0:
        assignments_df = assignments_df.join(
            membership_flags.loc[shared_idx, flag_cols], how='left'
        )
        assignments_df[flag_cols] = assignments_df[flag_cols].fillna(False)
        # Re-write with cohort flags
        assignments_df.to_csv(os.path.join(base_output, 'cluster_assignments.csv'))
        print(f"    Cluster assignments (with cohort flags) → cluster_assignments.csv")

    # (b) Cohort–cluster affinity
    affinity_df = compute_cohort_cluster_affinity(
        assignments_df, membership_flags, cohorts
    )
    affinity_path = os.path.join(base_output, 'cohort_cluster_affinity.csv')
    affinity_df.to_csv(affinity_path)
    print(f"    Cohort–cluster affinity → cohort_cluster_affinity.csv")
    print(affinity_df.to_string())

    # (c) Weighted feature list
    wf_df = compute_weighted_feature_list(H_df, feature_metadata)
    wf_df.to_csv(os.path.join(base_output, 'weighted_feature_list.csv'), index=False)
    print(f"    Weighted feature list   → weighted_feature_list.csv")

    # (d) Dominant cohort per cluster
    dom_df = annotate_cluster_dominant_cohort(affinity_df)
    dom_df.to_csv(os.path.join(base_output, 'cluster_dominant_cohort.csv'), index=False)
    print(f"    Dominant cohort map     → cluster_dominant_cohort.csv")
    for _, row in dom_df.iterrows():
        print(f"      {row['cluster']:20s} → {row['dominant_cohort']}  "
              f"({row['proportion']:.1%})")

    # (e) Affinity heatmap
    plot_cohort_cluster_affinity(
        affinity_df,
        os.path.join(base_fig, 'cohort_cluster_affinity.png'),
    )

    print("\n" + "=" * 60)
    print("COMBINED bNMF COMPLETE")
    print("=" * 60)
    print(f"  Optimal k  : {result['optimal_k']}")
    print(f"  Individuals: {result['n']}")
    print(f"  Features   : {H_df.shape[1]}")
    print(f"  Outputs  → {base_output}/")
    print(f"  Figures  → {base_fig}/")

    return result


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "bNMF subtype analysis per cohort using high-risk cases + "
            "low-control controls and features from cohort analysis outputs."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # ── Required ────────────────────────────────────────────────────────────
    parser.add_argument(
        "--pheno_data",
        help="Path to phenotype results directory (or PHENO_DATA env var)",
    )
    parser.add_argument(
        "--raw_features_file",
        help="Path to raw (untransformed) clinical features CSV "
             "(or RAW_FEATURES_FILE env var)",
    )

    # ── Feature selection ──────────────────────────────────────────────────
    parser.add_argument(
        "--filter_strategy", default='tiered',
        choices=['strict', 'majority', 'differential', 'tiered'],
        help="Genomic feature specificity filter (default: tiered)",
    )
    parser.add_argument(
        "--use_set", default='validation',
        choices=['holdout', 'validation', 'all'],
        help="Individual set used for bNMF discovery (default: validation); "
             "holdout is preserved for replication",
    )
    parser.add_argument(
        "--min_effect_size", type=float, default=0.0,
        help="Min |rank-biserial r| for clinical features (default: 0.0 = all significant)",
    )
    parser.add_argument(
        "--specificity_tier_max", type=int, default=3,
        choices=[1, 2, 3],
        help="Max genomic specificity tier to include (default: 3 = all tiers)",
    )

    # ── Individual selection ───────────────────────────────────────────────
    parser.add_argument(
        "--high_risk_threshold", type=int, default=HIGH_RISK_BIN_THRESHOLD,
        help=(
            f"bin > this value = high-risk case (default: {HIGH_RISK_BIN_THRESHOLD} "
            f"= top 20%%, deciles 9-10). Applied to cohort-specific bin column."
        ),
    )
    parser.add_argument(
        "--low_control_threshold", type=int, default=LOW_CONTROL_BIN_THRESHOLD,
        help=(
            f"bin < this value = low-control (default: {LOW_CONTROL_BIN_THRESHOLD} "
            f"= bottom 20%%, deciles 1-2). Must match the threshold used during "
            f"cohort clinical/genomic feature selection so the same control "
            f"population is used throughout."
        ),
    )

    # ── NMF parameters ─────────────────────────────────────────────────────
    parser.add_argument(
        "--k_min", type=int, default=2,
        help="Minimum k for NMF consensus (default: 2)",
    )
    parser.add_argument(
        "--k_max", type=int, default=6,
        help="Maximum k for NMF consensus (default: 6)",
    )
    parser.add_argument(
        "--n_runs", type=int, default=30,
        help="NMF restarts per k for consensus clustering (default: 30)",
    )
    parser.add_argument(
        "--k_select_method", default='elbow',
        choices=['elbow', 'max'],
        help="Optimal k selection method (default: elbow)",
    )
    parser.add_argument(
        "--l1_ratio", type=float, default=0.1,
        help="L1 ratio for NMF sparsity regularisation (default: 0.1)",
    )
    parser.add_argument(
        "--alpha_W", type=float, default=0.1,
        help="Alpha_W regularisation on W (basis matrix) (default: 0.1)",
    )
    parser.add_argument(
        "--alpha_H", type=float, default=0.1,
        help=(
            "Alpha_H regularisation on H (feature loadings) (default: 0.1). "
            "Must be > 0 when features include sparse/binary-like columns (e.g. "
            "HLA allele contributions). Setting alpha_H=0 allows H values to grow "
            "unconstrained, causing cluster collapse where everyone is assigned to "
            "cluster_1."
        ),
    )

    # ── Cohort selection ───────────────────────────────────────────────────
    parser.add_argument(
        "--cohorts", nargs='+', default=None,
        help="Specific cohort names to run (default: all found in PRS bins)",
    )
    parser.add_argument(
        "--no_holdout_prs", action='store_true',
        help="Exclude holdout PRS files when loading weighted feature values",
    )

    # ── Analysis mode ──────────────────────────────────────────────────────
    parser.add_argument(
        "--mode", default='both',
        choices=['per_cohort', 'combined', 'both', 'pca'],
        help=(
            "Analysis mode (default: both):\n"
            "  per_cohort — run bNMF separately within each cohort\n"
            "  combined   — run a single bNMF on the union of all cohorts\n"
            "  both       — run per-cohort then combined\n"
            "  pca        — build combined feature matrix and run PCA only;\n"
            "               no NMF. Assesses intrinsic dimensionality before\n"
            "               committing to NMF. Outputs scree plot, PC scatter\n"
            "               coloured by cohort, and top feature loadings per PC."
        ),
    )
    parser.add_argument(
        "--n_pca_components", type=int, default=30,
        help="Number of PCA components to compute in --mode pca (default: 30)",
    )
    parser.add_argument(
        "--clinical_only", action='store_true',
        help="PCA mode only: restrict matrix to clin_* columns (no genomic features). "
             "Useful for assessing clinical feature dimensionality independently.",
    )
    parser.add_argument(
        "--population", default='both',
        choices=['both', 'high_risk', 'low_control', 'separate'],
        help=(
            "Which population to use for per-cohort bNMF (default: both):\n"
            "  both        — union of high-risk cases and low-control controls\n"
            "  high_risk   — only high-risk cases (bin > threshold, PHENOTYPE=case)\n"
            "  low_control — only low-risk controls (bin < threshold, PHENOTYPE=ctrl)\n"
            "  separate    — run two independent bNMF analyses per cohort,\n"
            "                saving to {cohort}/high_risk/ and {cohort}/low_control/\n"
            "                so each population gets its own clusters and enrichment"
        ),
    )

    # ── Combined-mode feature options ─────────────────────────────────────
    parser.add_argument(
        "--weight_features", action='store_true',
        help="Enable importance weighting of features after scaling/filtering "
             "(weights features by effect_size_r / matching_z_score; default: off)",
    )
    parser.add_argument(
        "--include_prs_with_genomic", action='store_true',
        help=(
            "Include the cohort PRS anchor column even for cohorts that already "
            "have gen_ (SNP-level) features. Not recommended: PRS ≈ Σ(gen_ contributions) "
            "so including both introduces near-collinearity that can cause NMF cluster "
            "collapse. Default: PRS anchor is suppressed when gen_ features are present "
            "for that cohort; only added for cohorts with no selected SNP features."
        ),
    )
    parser.add_argument(
        "--scale_method", choices=['quantile', 'minmax'], default='quantile',
        help="Feature scaling before NMF. 'quantile' (default) maps each feature "
             "to uniform [0,1] by rank — removes the shared mean direction that "
             "causes collapse in homogeneous subgroups. 'minmax' preserves raw "
             "value distributions but amplifies correlated features.",
    )
    parser.add_argument(
        "--no_raw_clinical", action='store_true',
        help="Do not add raw clinical features to the combined matrix. "
             "By default all numeric clinical features are included alongside "
             "genomic features to provide real per-individual signal that is "
             "not subject to combine_first imputation artifacts.",
    )
    parser.add_argument(
        "--min_variance", type=float, default=0.005,
        help="Drop features with variance < this after scaling (default: 0.005)",
    )
    parser.add_argument(
        "--max_zero_fraction", type=float, default=0.80,
        help="Drop features where > this fraction of values are zero (default: 0.80)",
    )
    parser.add_argument(
        "--max_features", type=int, default=500,
        help="Cap total features by variance before NMF (default: 500)",
    )

    parser.add_argument(
        "--figures_only", action='store_true',
        help=(
            "Skip NMF entirely — reload saved cluster_profile.csv and "
            "feature_loadings.csv from existing result directories and "
            "regenerate only the PNG figures. Useful for updating plot "
            "styles without re-running the full analysis. "
            "Only --pheno_data is required in this mode."
        ),
    )
    parser.add_argument(
        "--top_n_features", type=int, default=20,
        help="Number of top features to show in cluster profile and loadings plots (default: 20)",
    )

    args = parser.parse_args()

    pheno_data        = args.pheno_data        or os.environ.get("PHENO_DATA")
    raw_features_file = args.raw_features_file or os.environ.get("RAW_FEATURES_FILE")

    if not pheno_data:
        raise ValueError("Provide --pheno_data or set the PHENO_DATA env var")

    # ── Figures-only shortcut ────────────────────────────────────────────────
    if args.figures_only:
        print("=" * 60)
        print("REPLOT FIGURES FROM SAVED CSVs")
        print("=" * 60)
        replot_figures_from_saved_csvs(pheno_data, top_n=args.top_n_features)
        raise SystemExit(0)

    if not raw_features_file:
        raise ValueError("Provide --raw_features_file or set the RAW_FEATURES_FILE env var")

    shared_kwargs = dict(
        pheno_data=pheno_data,
        raw_features_file=raw_features_file,
        filter_strategy=args.filter_strategy,
        use_set=args.use_set,
        k_min=args.k_min,
        k_max=args.k_max,
        n_runs=args.n_runs,
        k_select_method=args.k_select_method,
        l1_ratio=args.l1_ratio,
        alpha_W=args.alpha_W,
        alpha_H=args.alpha_H,
        min_effect_size=args.min_effect_size,
        specificity_tier_max=args.specificity_tier_max,
        high_risk_threshold=args.high_risk_threshold,
        low_control_threshold=args.low_control_threshold,
        cohorts_to_run=args.cohorts,
        include_holdout_prs=not args.no_holdout_prs,
        min_variance=args.min_variance,
        max_zero_fraction=args.max_zero_fraction,
        max_features=args.max_features,
        include_raw_clinical=not args.no_raw_clinical,
        scale_method=args.scale_method,
    )

    if args.mode == 'pca':
        run_pca_analysis(
            pheno_data=pheno_data,
            raw_features_file=raw_features_file,
            filter_strategy=args.filter_strategy,
            use_set=args.use_set,
            min_effect_size=args.min_effect_size,
            specificity_tier_max=args.specificity_tier_max,
            high_risk_threshold=args.high_risk_threshold,
            low_control_threshold=args.low_control_threshold,
            cohorts_to_run=args.cohorts,
            include_holdout_prs=not args.no_holdout_prs,
            min_variance=args.min_variance,
            max_zero_fraction=args.max_zero_fraction,
            max_features=args.max_features,
            include_prs_with_genomic=args.include_prs_with_genomic,
            population=args.population,
            include_raw_clinical=not args.no_raw_clinical,
            scale_method=args.scale_method,
            n_pca_components=args.n_pca_components,
            clinical_only=args.clinical_only,
        )
    else:
        if args.mode in ('per_cohort', 'both'):
            run_all_cohorts(**shared_kwargs, population=args.population)

        if args.mode in ('combined', 'both'):
            run_combined_bnmf(
                **shared_kwargs,
                weight_features=args.weight_features,
                include_prs_with_genomic=args.include_prs_with_genomic,
                population=args.population,
                clinical_only=args.clinical_only,
            )
