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

    print(f"    Individuals — high-risk: {stats['n_high_risk']:,}  "
          f"low-control: {stats['n_low_control']:,}  "
          f"total: {stats['n_total']:,}")
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
    Keep only PRS model entries whose name corresponds to a valid cohort.
    Drops covariates, all_*, all+env_combined, and any other model that does
    not match a cohort name — mirroring the use_all filter in the clinical
    comparison pipeline.
    """
    filtered = {
        k: v for k, v in prs_dict.items()
        if any(
            c.lower() == k.lower() or
            c.lower() in k.lower() or
            k.lower() in c.lower()
            for c in cohorts
        )
    }
    excluded = sorted(set(prs_dict.keys()) - set(filtered.keys()))
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
):
    """
    Build the feature matrix for bNMF for a single cohort.

    Combines:
      Clinical features — significant for this cohort (raw values from
                          raw_clinical_df), prefixed 'clin_'
      Genomic features  — SHAP-important for this cohort (per-individual
                          weighted contributions from prs_dict[model]),
                          prefixed 'gen_'

    Parameters
    ----------
    cohort                : str
    selected_iids         : pd.Index of IIDs to include
    genomic_features_df   : output of load_cohort_genomic_features()
    clinical_features_df  : output of load_cohort_clinical_features()
    prs_dict              : output of load_all_prs_feature_files()
    raw_clinical_df       : output of load_raw_clinical_data()
    min_effect_size       : float, minimum |effect_size_r| filter for clinical
    specificity_tier_max  : int, maximum genomic specificity tier to include

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
    # 2. Genomic weighted-contribution features
    # ------------------------------------------------------------------
    if not genomic_features_df.empty and 'cohort' in genomic_features_df.columns:
        cohort_gen = genomic_features_df[
            genomic_features_df['cohort'] == cohort
        ].copy()

        if 'specificity_tier' in cohort_gen.columns:
            cohort_gen = cohort_gen[
                cohort_gen['specificity_tier'] <= specificity_tier_max
            ]

        gen_features = cohort_gen['feature'].unique().tolist()

        if gen_features:
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
                available_gen = [f for f in gen_features if f in prs_feature_cols]

                if available_gen:
                    gen_sub = prs_df.loc[
                        prs_df.index.isin(selected_iids), available_gen
                    ].copy()
                    gen_sub.columns = [f'gen_{c}' for c in gen_sub.columns]
                    frames.append(gen_sub)
                    feat_counts['n_genomic'] = len(available_gen)
                    n_missing = len(gen_features) - len(available_gen)
                    print(f"    Genomic   : {len(available_gen)} features used"
                          f" (model '{model_key}')"
                          + (f" ({n_missing} not in PRS columns)" if n_missing else ""))
                else:
                    print(f"    [WARN] No genomic features matched PRS columns "
                          f"for cohort '{cohort}'")

    # ------------------------------------------------------------------
    # 3. Merge all feature blocks
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
    """Line plot of cophenetic correlation vs k."""
    if coph_df.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(coph_df['k'], coph_df['cophenetic_correlation'],
            marker='o', color='#1f77b4', linewidth=2, markersize=6)
    ax.set_xlabel('Number of clusters (k)')
    ax.set_ylabel('Cophenetic correlation')
    ax.set_title(f'Consensus NMF Stability — {cohort}')
    ax.set_xticks(coph_df['k'])
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
    Heatmap of per-cluster median profiles (top N features by variance
    across clusters).
    """
    if profile_df.empty:
        return
    median_cols = [c for c in profile_df.columns
                   if c.endswith('_median') and c not in ('n_median',)]
    if not median_cols:
        return

    # characterize_clusters() already uses 'cluster' as index; handle both cases
    if profile_df.index.name == 'cluster':
        profile_matrix = profile_df[median_cols].copy()
    elif 'cluster' in profile_df.columns:
        profile_matrix = profile_df.set_index('cluster')[median_cols].copy()
    else:
        profile_matrix = profile_df[median_cols].copy()
    profile_matrix.columns = [c.replace('_median', '') for c in median_cols]

    if len(profile_matrix.columns) > top_n:
        top_feats = profile_matrix.var(axis=0).nlargest(top_n).index.tolist()
        profile_matrix = profile_matrix[top_feats]

    fig_h = max(4, len(profile_matrix.columns) * 0.35)
    fig_w = max(4, len(profile_matrix) * 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        profile_matrix.T, cmap='RdBu_r', center=0,
        ax=ax, annot=True, fmt='.2f', annot_kws={'size': 7},
        linewidths=0.5, cbar_kws={'label': 'Median (scaled)'},
    )
    ax.set_title(f'Cluster Profiles (top {len(profile_matrix.columns)} features) — {cohort}')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Feature')
    plt.tight_layout()
    _save_fig(fig, output_path)


def plot_feature_loadings(H_df, output_path, cohort, top_n=20):
    """
    Per-cluster horizontal bar charts of top feature loadings.
    Clinical features are shown in red, genomic in blue.
    """
    if H_df.empty:
        return

    n_clusters = len(H_df)
    fig, axes = plt.subplots(
        1, n_clusters,
        figsize=(n_clusters * 5, max(4, min(top_n * 0.35, 12))),
    )
    if n_clusters == 1:
        axes = [axes]

    for ax, (cluster_label, row) in zip(axes, H_df.iterrows()):
        top = row.nlargest(top_n)
        colours = ['#d62728' if f.startswith('clin_') else '#1f77b4'
                   for f in top.index]
        labels  = [f.replace('clin_', 'C: ').replace('gen_', 'G: ')
                   for f in top.index]

        ax.barh(range(len(top)), top.values, color=colours, alpha=0.85)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_title(cluster_label, fontsize=9)
        ax.set_xlabel('Loading')

    legend_handles = [
        mpatches.Patch(facecolor='#d62728', label='Clinical'),
        mpatches.Patch(facecolor='#1f77b4', label='Genomic'),
    ]
    axes[-1].legend(handles=legend_handles, loc='lower right', fontsize=8)
    plt.suptitle(f'Top Feature Loadings — {cohort}', fontsize=11)
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

    Individuals are the UNION of high-risk cases and low-control controls across all cohorts.

    Returns
    -------
    combined_df      : pd.DataFrame  (union_IIDs × deduplicated features)
    feature_metadata : dict  {column_name: {'cohort_origin', 'type', 'weight'}}
    per_cohort_iids  : dict  {cohort: pd.Index of selected IIDs}
    """
    per_cohort_iids = {}
    all_iids        = pd.Index([])

    # Collect per-cohort individual sets and union of all IIDs
    for cohort in cohorts:
        iids, _ = select_cohort_individuals(
            prs_bins, cohort,
            high_risk_threshold=high_risk_threshold,
            low_control_threshold=low_control_threshold,
        )
        per_cohort_iids[cohort] = iids
        all_iids = all_iids.union(iids)

    if all_iids.empty:
        return pd.DataFrame(), {}, per_cohort_iids

    feature_metadata = {}
    frames           = []

    # ── 1. Genomic features — collected first so we know which cohorts they
    #       cover before deciding which PRS anchors are needed.
    #       Deduplicated: one column per SNP, weight = max matching_z_score.
    # ──────────────────────────────────────────────────────────────────────
    gen_best: dict[str, dict] = {}   # feature_name → {z_score, cohort, model_key}
    if not genomic_features_df.empty and 'cohort' in genomic_features_df.columns:
        for cohort in cohorts:
            g_sub = genomic_features_df[
                genomic_features_df['cohort'] == cohort
            ].copy()
            if 'specificity_tier' in g_sub.columns:
                g_sub = g_sub[g_sub['specificity_tier'] <= specificity_tier_max]
            model_key = _find_model_key(cohort, prs_dict)
            if model_key is None:
                continue
            z_col = 'matching_z_score' if 'matching_z_score' in g_sub.columns else None
            for f in g_sub['feature'].unique():
                z = abs(float(g_sub.loc[g_sub['feature'] == f, z_col].iloc[0])) \
                    if z_col else 0.0
                if f not in gen_best or z > gen_best[f]['z_score']:
                    gen_best[f] = {
                        'z_score':       z,
                        'cohort_origin': cohort,
                        'model_key':     model_key,
                    }

    # Track which cohorts are represented by at least one gen_ feature.
    # For those cohorts the PRS anchor is redundant (PRS ≈ Σ gen_ contributions)
    # and is suppressed unless include_prs_with_genomic=True.
    cohorts_with_genomic = {meta['cohort_origin'] for meta in gen_best.values()}

    if gen_best:
        from collections import defaultdict
        model_features: dict = defaultdict(list)
        for f, meta in gen_best.items():
            model_features[meta['model_key']].append(f)

        gen_frames = []
        for model_key, feats in model_features.items():
            prs_df = prs_dict[model_key]
            valid  = [f for f in feats
                      if f in prs_df.columns
                      and f not in PRS_META_COLS
                      and not f.startswith('PC')]
            if not valid:
                continue
            blk = prs_df.loc[prs_df.index.isin(all_iids), valid].copy()
            blk.columns = [f'gen_{f}' for f in blk.columns]
            gen_frames.append(blk)

        if gen_frames:
            gen_block = gen_frames[0]
            for extra in gen_frames[1:]:
                gen_block = gen_block.join(extra, how='outer')
            frames.append(gen_block)
            for f, meta in gen_best.items():
                col = f'gen_{f}'
                if col in gen_block.columns:
                    feature_metadata[col] = {
                        'cohort_origin': meta['cohort_origin'],
                        'type':          'genomic',
                        'weight':        meta['z_score'],
                    }
            print(f"    Genomic   : {len(gen_best)} deduplicated features "
                  f"(best model per SNP by matching_z_score)")
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
    max_zero_fraction=0.90,
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
    output_path=None,
    fig_path=None,
    min_variance=0.005,
    max_zero_fraction=0.90,
    max_features=500,
):
    """
    Run consensus bNMF for one cohort and write all outputs.

    Parameters
    ----------
    min_variance     : float — drop columns with variance below this after scaling
    max_zero_fraction: float — drop columns that are > this fraction zero
    max_features     : int   — cap total features by variance (combined mode only)

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
    X_df, kept_cols = bnmf_core.impute_and_scale(feature_matrix)

    # Remove near-constant and highly sparse columns before NMF
    X_df, kept_cols = _variance_sparsity_filter(
        X_df,
        min_variance=min_variance,
        max_zero_fraction=max_zero_fraction,
        max_features=max_features,
    )
    X = X_df.values

    if X.shape[1] < 2:
        print(f"  [SKIP] {cohort}: fewer than 2 features after preprocessing")
        return None

    print(f"\n  Consensus bNMF: {X.shape[0]:,} individuals × {X.shape[1]} features")

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
        plot_consensus_map(
            k_results[optimal_k]['C'], optimal_k,
            os.path.join(fig_path, f'consensus_map_k{optimal_k}.png'),
            cohort,
        )
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
# MAIN ORCHESTRATOR
# ============================================================================

def run_all_cohorts(
    pheno_data,
    raw_features_file,
    filter_strategy='tiered',
    use_set='holdout',
    k_min=2,
    k_max=6,
    n_runs=30,
    k_select_method='elbow',
    l1_ratio=0.1,
    alpha_W=0.1,
    min_effect_size=0.0,
    specificity_tier_max=3,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    cohorts_to_run=None,
    include_holdout_prs=True,
    population='both',
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
                output_path=cohort_output,
                fig_path=cohort_figs,
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
# COMBINED ORCHESTRATOR
# ============================================================================

def run_combined_bnmf(
    pheno_data,
    raw_features_file,
    filter_strategy='tiered',
    use_set='holdout',
    k_min=2,
    k_max=6,
    n_runs=30,
    k_select_method='elbow',
    l1_ratio=0.1,
    alpha_W=0.1,
    min_effect_size=0.0,
    specificity_tier_max=3,
    high_risk_threshold=HIGH_RISK_BIN_THRESHOLD,
    low_control_threshold=LOW_CONTROL_BIN_THRESHOLD,
    cohorts_to_run=None,
    include_holdout_prs=True,
    weight_features=True,
    min_variance=0.005,
    max_zero_fraction=0.90,
    max_features=500,
    include_prs_with_genomic=False,
):
    """
    Combined cross-cohort bNMF subtype analysis.

    Builds a single de-duplicated feature matrix:
      {cohort}_prs   — scaled_prs per cohort model (non-redundant anchors)
      clin_{feature} — raw clinical value, one column; weight = max |effect_size_r|
      gen_{feature}  — weighted contribution from best model; weight = max matching_z_score

    This avoids the cluster-collapse caused by cohort-prefixed redundant columns.

    Individuals are the UNION of high-risk cases and low-control controls across all cohorts.

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
    base_output = os.path.join(scores_path,  'combinedCohortBnmf')
    base_fig    = os.path.join(pheno_data,   'figures', 'combinedCohortBnmf')

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

    # ── Build combined feature matrix (de-duplicated) ─────────────────────
    print("\n[4] Building de-duplicated feature matrix...")
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

    # ── Apply feature importance weighting ────────────────────────────────
    # Multiply each column by its importance weight before imputation+scaling
    # so the NMF sees stronger signal in high-importance features.
    # PRS anchors (weight=1.0) are left unchanged.
    nmf_matrix = combined_matrix.copy()
    if weight_features and feature_metadata:
        for col in nmf_matrix.columns:
            w = feature_metadata.get(col, {}).get('weight', 1.0)
            if w > 0:
                nmf_matrix[col] = nmf_matrix[col] * w
        print(f"    Feature importance weighting applied "
              f"({sum(1 for m in feature_metadata.values() if m['weight'] > 0)} "
              f"features with weight > 0)")

        # Save weights alongside outputs
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
        "--use_set", default='holdout',
        choices=['holdout', 'validation', 'all'],
        help="Individual set used for clinical comparison output (default: holdout)",
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
        help="Alpha_W regularisation parameter (default: 0.1)",
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
        choices=['per_cohort', 'combined', 'both'],
        help=(
            "Analysis mode (default: both):\n"
            "  per_cohort — run bNMF separately within each cohort\n"
            "  combined   — run a single bNMF on the union of all cohorts,\n"
            "               producing comparable clusters + pathway outputs\n"
            "  both       — run per-cohort then combined"
        ),
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
        "--no_weight_features", action='store_true',
        help="Disable importance weighting of features before NMF "
             "(default: features are weighted by effect_size_r / matching_z_score)",
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
        "--min_variance", type=float, default=0.005,
        help="Drop features with variance < this after scaling (default: 0.005)",
    )
    parser.add_argument(
        "--max_zero_fraction", type=float, default=0.90,
        help="Drop features where > this fraction of values are zero (default: 0.90)",
    )
    parser.add_argument(
        "--max_features", type=int, default=500,
        help="Cap total features by variance before NMF (default: 500)",
    )

    args = parser.parse_args()

    pheno_data        = args.pheno_data        or os.environ.get("PHENO_DATA")
    raw_features_file = args.raw_features_file or os.environ.get("RAW_FEATURES_FILE")

    if not pheno_data:
        raise ValueError("Provide --pheno_data or set the PHENO_DATA env var")
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
        min_effect_size=args.min_effect_size,
        specificity_tier_max=args.specificity_tier_max,
        high_risk_threshold=args.high_risk_threshold,
        low_control_threshold=args.low_control_threshold,
        cohorts_to_run=args.cohorts,
        include_holdout_prs=not args.no_holdout_prs,
    )

    if args.mode in ('per_cohort', 'both'):
        run_all_cohorts(**shared_kwargs, population=args.population)

    if args.mode in ('combined', 'both'):
        run_combined_bnmf(
            **shared_kwargs,
            weight_features=not args.no_weight_features,
            min_variance=args.min_variance,
            max_zero_fraction=args.max_zero_fraction,
            max_features=args.max_features,
            include_prs_with_genomic=args.include_prs_with_genomic,
        )
