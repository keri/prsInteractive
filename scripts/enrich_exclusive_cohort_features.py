#!/usr/bin/env python3
"""
enrich_exclusive_cohort_features.py

Enriches important features from high-risk cases EXCLUSIVE to each statistically
distinct cohort against functional cluster annotations.

Pipeline:
  1. Identify exclusive high-risk cases per cohort: individuals in the top 20%
     of the risk bin in ONE cohort only (not high-risk in any other), PHENOTYPE == 2.
     Controls are within-cohort bottom-20% low-risk individuals (PHENOTYPE == 1).
     Bin thresholds are derived proportionally from the scale in the loaded file
     (1–10 deciles for training; 1–1000 quantiles for holdout) so the same
     top-20% / bottom-20% split is preserved in both cases.

  2. Find important features discriminating exclusive cases from controls:
       Clinical features : Mann-Whitney U on raw_features_file
                           (unbalanced = all controls; balanced = bootstrap n_cases
                            controls over --n_iterations iterations)
       Genomic features  : loaded from pre-computed SHAP output
                           (importantFeaturesAcrossCohortsAndTrainingData.Filtered.*.csv)
                           — already cohort-exclusive by cross-cohort filtering strategy

  3. Merge genomic features with feature_scores_annotated.csv to get snp1/snp2
     and their gene names. GxG features (two rsIDs) expand to one row per SNP.

  4. Load functional cluster CSVs from {results_path}/subtypeResearch/:
       - pattern: *functional_cluster*.csv
       - expected feature identifier cols: feature, Feature, Trait, trait, Loci, loci
       - if "cluster assignment" / "cluster_assignment" column is present, pivot
         long → wide using beta / effect as values

  5. Map each feature to functional clusters (binary 0/1):
       Priority 1 – feature → feature : SNP rsID or clinical name vs cluster feature col
       Priority 2 – loci   → loci     : feature name (when it IS a loci token) vs loci col
       Priority 3 – feature → loci    : gene-name tokens from snp_gene_names vs loci col
       match_method column records which level produced the hit.

Output:
  {scores_path}/functionalEnrichment/
    functional_cluster_enrichment_{use_set}_{balance_strategy}.csv

    Columns:
      raw_feature        – original feature string from cohort analysis
      feature            – resolved SNP rsID or clinical measure (1 per row)
      cohort             – source cohort
      p_value            – Mann-Whitney p (clinical) or None (genomic SHAP)
      p_value_fdr        – BH-corrected p
      effect_size_r      – rank-biserial r (clinical) or odds_ratio (genomic)
      snp1_gene_names    – semicolon-joined gene names for snp1 (if annotated)
      snp2_gene_names    – semicolon-joined gene names for snp2 (if annotated)
      feature_source     – 'genomic_shap' or 'clinical'
      match_method       – 'feature_feature' | 'loci_loci' | 'feature_loci' | 'none'
      n_clusters_matched – count of functional cluster columns with value 1
      <cluster_col> ...  – one binary (0/1) column per functional cluster

Design note — balanced vs unbalanced:
  Unbalanced (default): use all controls. Maximises statistical power; FDR
  correction handles inflated discovery rate. Preferred when control pool is
  large (>> n_cases).
  Balanced: bootstrap n_cases controls per iteration, aggregate median p-value.
  Reduces influence of control sample composition; adds computational cost.
  Useful when control pool is only 2-5x cases, or when class imbalance drives
  spurious discoveries in validation.

Usage:
  python enrich_exclusive_cohort_features.py \\
    --pheno_data  /path/to/combinedAnalysis \\
    --raw_features_file /path/to/raw_features.csv \\
    [--filter_strategy tiered] \\
    [--use_set holdout] \\
    [--balance_strategy unbalanced] \\
    [--n_iterations 1000] \\
    [--alpha 0.05] \\
    [--min_effect_size 0.1] \\
    [--use_all]
"""

import pandas as pd
import numpy as np
import os
import glob
import argparse
import re
import warnings

from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# =============================================================================
#  COLUMN-NAME NORMALISATION
#  Collapses variant spellings and separators so that names like
#  "residual glycaemia", "residual_glycaemic", "residual glycemia"
#  all resolve to the same canonical key before any matching is done.
# =============================================================================

# Runtime-populated clinical feature map: {normalised_alias → canonical_key}.
# Loaded from clinical_feature_map.csv in subtypeResearch/.  Empty until
# load_clinical_feature_map() is called in main().
_CLINICAL_FEATURE_MAP: dict = {}

# (pattern, replacement) — applied left-to-right after basic cleaning.
# American / informal → canonical British / clinical spelling used in UKB.
_SPELLING_CANON = [
    # glyc*
    (r'glycemi',   'glycaemi'),     # glycemia → glycaemia, glycemic → glycaemic
    (r'glycaemic', 'glycaemia'),    # adjective → noun stem
    # haem / hem
    (r'hemoglob',  'haemoglob'),
    (r'hematocr',  'haematocr'),
    (r'hematin',   'haematin'),
    # anaemia / anemia
    (r'\baneми',   'anaemi'),       # unicode-safe stub (real below)
    (r'\banemi',   'anaemi'),
    # oestr* / estr*
    (r'\bestrog',  'oestrog'),
    # normalise common trailing noise: trailing 'ic' adjective → stem
    # e.g. glycaemic → glycaemia already handled above
    # cholesterol variants
    (r'ldl_c\b',   'ldl_cholesterol'),
    (r'hdl_c\b',   'hdl_cholesterol'),
    # fasting variants
    (r'fasting_glucose',   'fasting_glucose'),   # keep canonical
    (r'fast_glucose',      'fasting_glucose'),
    (r'fasting_blood_glucose', 'fasting_glucose'),
    # hba1c / a1c / glycated haemoglobin
    (r'\bhba1c\b', 'hba1c'),
    (r'\bhaemoglobin_a1c\b', 'hba1c'),
    (r'\bhemoglobin_a1c\b',  'hba1c'),
    # BMI variants
    (r'body_mass_index', 'bmi'),
    # residual b-cell
    (r'residual_beta[_-]?cell', 'residual_beta_cell'),
]


def normalise_name(name: str) -> str:
    """
    Produce a canonical string key for a feature / column name.

    Steps
    -----
    1. Lowercase and strip surrounding whitespace.
    2. Replace any run of whitespace, hyphens, forward-slashes or dots with '_'.
    3. Collapse repeated underscores; strip leading/trailing underscores.
    4. Apply known spelling-equivalence substitutions in _SPELLING_CANON.

    Examples
    --------
    "Residual Glycaemia"  → "residual_glycaemia"
    "residual_glycaemic"  → "residual_glycaemia"   (adjective → noun stem)
    "residual glycemia"   → "residual_glycaemia"   (US spelling)
    "HbA1c"               → "hba1c"
    "LDL-C"               → "ldl_cholesterol"
    """
    s = _normalise_base(name)
    if not s:
        return ''
    # Clinical feature map lookup (populated at runtime from clinical_feature_map.csv).
    # Resolves aliases like 'hba1c', 'glycated_hemoglobin' → canonical primary term.
    if _CLINICAL_FEATURE_MAP:
        s = _CLINICAL_FEATURE_MAP.get(s, s)
    return s


def normalise_series_names(col_names):
    """
    Return a dict mapping each original column name to its canonical form,
    and a reverse dict mapping canonical → list of originals (for diagnostics).
    """
    canon_map = {c: normalise_name(c) for c in col_names}
    # Warn about collisions (different raw names → same canonical)
    reverse: dict = {}
    for orig, can in canon_map.items():
        reverse.setdefault(can, []).append(orig)
    collisions = {k: v for k, v in reverse.items() if len(v) > 1}
    if collisions:
        for canon, originals in collisions.items():
            print(f"  [NORM] Collapsed to '{canon}': {originals}")
    return canon_map, reverse


def load_clinical_feature_map(subtype_path: str) -> dict:
    """
    Load ``clinical_feature_map.csv`` from *subtype_path* and build a
    normalised-alias → canonical-key lookup.

    CSV structure
    -------------
    feature   – primary / canonical name (space-separated; most common form
                used in the functional cluster files).  This becomes the
                canonical key after normalisation.
    feature1  – first alternative spelling / abbreviation.
    feature2  – second alternative spelling / abbreviation.

    Any number of extra columns are silently ignored.

    Logic
    -----
    1. Normalise the *feature* value (lowercase, separator→_, spelling canon)
       → this is the canonical key.
    2. For each non-empty cell in feature1 / feature2, normalise that value
       and store alias → canonical_key in the returned dict.
    3. The canonical key itself is also stored (identity mapping) so that
       names that are already in canonical form pass through unchanged.

    Returns
    -------
    dict  {normalised_alias: normalised_canonical}

    The returned dict is merged into the module-level ``_CLINICAL_FEATURE_MAP``
    so that ``normalise_name()`` can consult it on every call.
    """
    map_path = os.path.join(subtype_path, 'clinical_feature_map.csv')
    if not os.path.isfile(map_path):
        print(f"  [WARN] clinical_feature_map.csv not found at {map_path}; "
              "clinical alias expansion disabled.")
        return {}

    df = pd.read_csv(map_path, encoding='utf-8-sig')
    # Drop completely empty columns (the CSV has trailing empty cols)
    df = df.dropna(axis=1, how='all')
    df.columns = df.columns.str.strip()

    alias_map: dict = {}
    n_entries = 0

    for _, row in df.iterrows():
        raw_primary = row.get('feature', '')
        if not raw_primary or (isinstance(raw_primary, float) and np.isnan(raw_primary)):
            continue
        # Normalise without consulting the (still-empty) map to avoid recursion
        canonical = _normalise_base(str(raw_primary))
        # Identity mapping: canonical → canonical
        alias_map[canonical] = canonical

        for alt_col in ('feature1', 'feature2'):
            raw_alt = row.get(alt_col, '')
            if not raw_alt or (isinstance(raw_alt, float) and np.isnan(raw_alt)):
                continue
            alt_norm = _normalise_base(str(raw_alt))
            if alt_norm and alt_norm != canonical:
                alias_map[alt_norm] = canonical
                n_entries += 1

    print(f"  Loaded clinical_feature_map.csv: {len(df)} entries, "
          f"{n_entries} alias mappings registered.")
    return alias_map


def _normalise_base(name: str) -> str:
    """
    Core normalisation without the clinical-map lookup (used when building
    the map itself to avoid circular reference).
    """
    if not name or (isinstance(name, float) and np.isnan(name)):
        return ''
    s = str(name).lower().strip()
    s = re.sub(r'[\s\-/\.]+', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    for pattern, repl in _SPELLING_CANON:
        s = re.sub(pattern, repl, s)
    return s


# =============================================================================
#  STEP 1: EXCLUSIVE HIGH-RISK CASE IDENTIFICATION
# =============================================================================

def detect_bin_thresholds(combined_prs, bin_cols):
    """
    Infer high-risk and low-risk bin thresholds from the actual scale in the
    loaded file, preserving the same 80th / 20th percentile split regardless
    of whether bins run 1–10 (deciles, training file) or 1–1000 (quantiles,
    holdout file).

    Equivalences
    ------------
    1–10  scale : high > 8  (top 20%),  low < 3   (bottom 20%)
    1–1000 scale: high > 800 (top 20%), low < 201  (bottom 20%)
    General     : high > 0.8 × max_bin,  low < 0.2 × max_bin + 1

    Returns
    -------
    (high_threshold, low_threshold) — both floats; compare with > and <
    """
    max_bin = float(combined_prs[bin_cols].max().max())
    high_threshold = 0.8 * max_bin
    low_threshold  = 0.2 * max_bin + 1
    print(f"  Bin scale detected: 1–{max_bin:.0f}  →  "
          f"high risk bin > {high_threshold:.0f}  |  "
          f"low risk bin < {low_threshold:.0f}")
    return high_threshold, low_threshold


def identify_exclusive_groups(combined_prs, bin_cols):
    """
    For each cohort identify:
      - exclusive high-risk cases : bin_{cohort} > high_threshold AND NOT above
                                     high_threshold in any other cohort, PHENOTYPE == 2
      - within-cohort controls    : bin_{cohort} < low_threshold AND PHENOTYPE == 1

    Thresholds are derived proportionally from the bin scale so that the same
    top-20% / bottom-20% split is applied whether the file uses 1–10 deciles
    (training set) or 1–1000 quantiles (holdout set).

    Parameters
    ----------
    combined_prs : DataFrame indexed by IID, columns include bin_* and PHENOTYPE
    bin_cols     : list of bin_* column names to consider

    Returns
    -------
    dict  {cohort: {'cases': [IID, ...], 'controls': [IID, ...]}}
    """
    high_thresh, low_thresh = detect_bin_thresholds(combined_prs, bin_cols)
    cohort_groups = {}
    cohort_names = [c.replace('bin_', '') for c in bin_cols]

    for col, cohort in zip(bin_cols, cohort_names):
        is_high_this = combined_prs[col] > high_thresh
        other_cols = [c for c in bin_cols if c != col]
        is_high_any_other = (
            combined_prs[other_cols].gt(high_thresh).any(axis=1)
            if other_cols
            else pd.Series(False, index=combined_prs.index)
        )

        exclusive_cases = combined_prs[
            is_high_this & ~is_high_any_other & (combined_prs['PHENOTYPE'] == 2)
        ].index.tolist()

        controls = combined_prs[
            (combined_prs[col] < low_thresh) & (combined_prs['PHENOTYPE'] == 1)
        ].index.tolist()

        cohort_groups[cohort] = {'cases': exclusive_cases, 'controls': controls}
        print(f"  {cohort}: {len(exclusive_cases):>5} exclusive high-risk cases  "
              f"| {len(controls):>6} controls")

    return cohort_groups


# =============================================================================
#  STEP 2a: CLINICAL FEATURE IMPORTANCE (Mann-Whitney U)
# =============================================================================

def _mwu(a_series, b_series):
    """Mann-Whitney U, return (p_value, rank_biserial_r) or (None, None)."""
    a = a_series.dropna().values
    b = b_series.dropna().values
    if len(a) < 5 or len(b) < 5:
        return None, None
    try:
        stat, pval = stats.mannwhitneyu(a, b, alternative='two-sided')
        r = 1 - (2 * stat) / (len(a) * len(b))
        return pval, r
    except Exception:
        return None, None


def _apply_fdr(result_df, p_col='p_value', alpha=0.05, min_effect=0.1):
    """Add p_value_fdr and significant columns to result_df in place."""
    valid = result_df[p_col].notna()
    pvals = result_df.loc[valid, p_col].values
    if len(pvals) == 0:
        result_df['p_value_fdr'] = np.nan
        result_df['significant'] = False
        return result_df
    _, fdr_pvals, _, _ = multipletests(pvals, method='fdr_bh')
    result_df.loc[valid, 'p_value_fdr'] = fdr_pvals
    result_df['p_value_fdr'] = result_df['p_value_fdr'].fillna(1.0)
    result_df['significant'] = (
        (result_df['p_value_fdr'] < alpha) &
        (result_df['effect_size_r'].abs() >= min_effect)
    )
    return result_df


def find_clinical_features_unbalanced(raw_df, cases, controls, feature_cols,
                                       alpha, min_effect_size):
    """
    Unbalanced: all controls vs exclusive cases.
    Statistical power is maximised; recommended when control pool >> n_cases.
    """
    df_cases = raw_df.loc[raw_df.index.isin(cases)]
    df_ctrl = raw_df.loc[raw_df.index.isin(controls)]

    rows = []
    for feat in feature_cols:
        pval, effect = _mwu(df_cases[feat], df_ctrl[feat])
        if pval is not None:
            rows.append({
                'feature': feat,
                'p_value': pval,
                'effect_size_r': effect,
                'n_cases': int(df_cases[feat].notna().sum()),
                'n_controls': int(df_ctrl[feat].notna().sum()),
            })

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    return _apply_fdr(result, alpha=alpha, min_effect=min_effect_size)


def find_clinical_features_balanced(raw_df, cases, controls, feature_cols,
                                     alpha, min_effect_size,
                                     n_iterations=1000, rng_seed=42):
    """
    Balanced: bootstrap n_cases controls per iteration, use median p-value.

    Feature significance: FDR-corrected median p < alpha AND |r| >= min_effect_size
    (effect size computed once on the full unbalanced comparison for stability).

    Use when control pool is only 2–5× the number of cases, or when imbalance
    is suspected to drive spurious clinical associations.
    """
    n_cases = len(cases)
    df_cases = raw_df.loc[raw_df.index.isin(cases)]
    df_ctrl_pool = raw_df.loc[raw_df.index.isin(controls)]

    replace = len(df_ctrl_pool) < n_cases
    if replace:
        print(f"    [WARN] Fewer controls ({len(df_ctrl_pool)}) than cases ({n_cases}); "
              "sampling with replacement.")

    rng = np.random.default_rng(rng_seed)
    pval_accumulator = {feat: [] for feat in feature_cols}

    for _ in range(n_iterations):
        idx = rng.choice(df_ctrl_pool.index, size=n_cases, replace=replace)
        df_sample = df_ctrl_pool.loc[idx]
        for feat in feature_cols:
            pval, _ = _mwu(df_cases[feat], df_sample[feat])
            if pval is not None:
                pval_accumulator[feat].append(pval)

    # Effect size from full (unbalanced) run for interpretability
    rows = []
    for feat in feature_cols:
        pvals = pval_accumulator[feat]
        if not pvals:
            continue
        _, effect = _mwu(df_cases[feat], df_ctrl_pool[feat])
        rows.append({
            'feature': feat,
            'p_value': float(np.median(pvals)),
            'p_value_median_across_iterations': float(np.median(pvals)),
            'effect_size_r': effect if effect is not None else np.nan,
            'n_cases': int(df_cases[feat].notna().sum()),
            'n_controls_pool': int(df_ctrl_pool[feat].notna().sum()),
            'n_iterations': len(pvals),
        })

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    return _apply_fdr(result, alpha=alpha, min_effect=min_effect_size)


# =============================================================================
#  STEP 2b: GENOMIC FEATURES (from pre-computed SHAP)
# =============================================================================

def load_genomic_features(scores_path, filter_strategy):
    """
    Load SHAP-based cohort-exclusive genomic features from pre-computed output.
    These are already filtered for cross-cohort exclusivity.

    Returns DataFrame with columns: feature, cohort, odds_ratio, [specificity_tier]
    """
    path = os.path.join(
        scores_path, 'importantCohortScores',
        f'importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv'
    )
    if not os.path.exists(path):
        print(f"  [WARN] Pre-computed SHAP features not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    # Normalise cohort column name
    if 'feature_set' in df.columns and 'cohort' not in df.columns:
        df.rename(columns={'feature_set': 'cohort'}, inplace=True)

    # Include both the strict cross-cohort RF-SHAP features (HighCases) and the
    # per-cohort exclusive bootstrap GBM-SHAP features (ExclusiveHighCases).
    # LowControls are excluded — the enrichment is case-focussed.
    if 'training_data_used' in df.columns:
        df = df[df['training_data_used'].isin(['HighCases', 'ExclusiveHighCases'])]

    print(f"  Loaded {len(df)} pre-computed genomic features "
          f"({df['cohort'].nunique()} cohorts)")
    if 'training_data_used' in df.columns:
        print('  ' + df.groupby('training_data_used').size()
              .rename('n').reset_index().to_string(index=False))
    return df


# =============================================================================
#  STEP 3: MERGE WITH feature_scores_annotated.csv AND EXPAND SNPs
# =============================================================================

def load_feature_annotation(scores_path):
    """Load feature_scores_annotated.csv if it exists."""
    path = os.path.join(scores_path, 'feature_scores_annotated.csv')
    if not os.path.exists(path):
        print(f"  [WARN] feature_scores_annotated.csv not found at: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    print(f"  Loaded feature_scores_annotated: {len(df)} rows")
    return df


def expand_to_snp_rows(feature_df, annotation_df):
    """
    Expand combined feature list to one row per SNP (or clinical measure).

    For G features   : one row, feature = snp1 (or raw feature if unannotated)
    For GxG features : two rows, one per snp (snp1 and snp2)
    For clinical     : one row, feature = raw_feature (no SNP)

    Parameters
    ----------
    feature_df   : DataFrame with columns raw_feature, cohort, feature_source,
                   and statistical cols (p_value, p_value_fdr, effect_size_r,
                   odds_ratio, specificity_tier — any subset)
    annotation_df: feature_scores_annotated DataFrame; can be empty

    Returns
    -------
    DataFrame with additional columns: feature, snp1, snp1_gene_names,
                                       snp2, snp2_gene_names
    """
    # Build fast lookup: raw feature string → annotation row
    ann_lookup = {}
    if not annotation_df.empty and 'feature' in annotation_df.columns:
        for _, arow in annotation_df.iterrows():
            ann_lookup[str(arow['feature'])] = arow

    expanded = []
    for _, row in feature_df.iterrows():
        raw = str(row['raw_feature'])
        base = row.to_dict()
        ann = ann_lookup.get(raw)

        snp1 = ann.get('snp1') if ann is not None and pd.notna(ann.get('snp1')) else None
        snp2 = ann.get('snp2') if ann is not None and pd.notna(ann.get('snp2')) else None
        snp1_genes = ann.get('snp1_gene_names') if ann is not None else None
        snp2_genes = ann.get('snp2_gene_names') if ann is not None else None

        # Row for snp1 (or raw feature when no annotation)
        resolved_feature = snp1 if snp1 else raw
        expanded.append({
            **base,
            'feature': resolved_feature,
            'snp1': snp1,
            'snp1_gene_names': snp1_genes,
            'snp2': snp2,
            'snp2_gene_names': snp2_genes,
        })

        # Additional row for snp2 in GxG features
        if snp2:
            expanded.append({
                **base,
                'feature': snp2,
                'snp1': snp1,
                'snp1_gene_names': snp1_genes,
                'snp2': snp2,
                'snp2_gene_names': snp2_genes,
            })

    return pd.DataFrame(expanded)


# =============================================================================
#  STEP 4: LOAD AND NORMALISE FUNCTIONAL CLUSTER FILES
# =============================================================================

def _first_col(df, candidates):
    """Return name of first candidate column found in df, else None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_cluster_file(filepath):
    """
    Load one functional-cluster CSV and return a normalised wide DataFrame.

    If the file has a 'cluster assignment' / 'cluster_assignment' column it is
    pivoted from long → wide (index = feature-identifier cols, columns = cluster
    assignment values, values = beta / effect).

    Returns
    -------
    wide_df      : DataFrame with identifier cols + one column per cluster
    cluster_cols : list of cluster column names in wide_df
    feat_col     : name of the primary feature/SNP identifier column (or None)
    loci_col     : name of the loci column (or None)
    trait_col    : name of the trait column (or None)
    """
    df = pd.read_csv(filepath)

    # 'feature' / 'Feature' covers both SNP rsID columns and Trait columns that
    # have already been renamed to 'feature' on the cluster-file side.
    feat_col  = _first_col(df, ['feature', 'Feature', 'SNP', 'snp', 'snp1'])
    # 'loci' / 'Loci' also covers snp1_gene_names / snp2_gene_names columns that
    # carry gene-name annotations equivalent to our snp1_gene_names field.
    loci_col  = _first_col(df, ['loci', 'Loci', 'Locus', 'locus',
                                 'snp1_gene_names', 'snp2_gene_names',
                                 'snp1_gene_name',  'snp2_gene_name'])
    trait_col = _first_col(df, ['Trait', 'trait'])
    ca_col    = _first_col(df, ['cluster assignment', 'cluster_assignment'])
    eff_col   = _first_col(df, ['beta', 'effect'])

    id_cols = [c for c in [feat_col, loci_col, trait_col] if c is not None]

    if ca_col is not None:
        # Long → wide
        if eff_col is None:
            df['_presence'] = 1
            eff_col = '_presence'

        pivot_index = id_cols if id_cols else list(df.columns[:1])
        try:
            wide = (
                df.pivot_table(
                    index=pivot_index,
                    columns=ca_col,
                    values=eff_col,
                    aggfunc='first',
                )
                .reset_index()
            )
            wide.columns = [str(c) for c in wide.columns]
            cluster_cols = [c for c in wide.columns if c not in id_cols]
            return wide, cluster_cols, feat_col, loci_col, trait_col
        except Exception as e:
            print(f"    [WARN] pivot failed for {os.path.basename(filepath)}: {e}; "
                  "treating whole file as single cluster")
            file_stem = os.path.splitext(os.path.basename(filepath))[0]
            df[file_stem] = 1
            return df[id_cols + [file_stem]], [file_stem], feat_col, loci_col, trait_col
    else:
        # No cluster assignment column — remaining non-identifier columns ARE
        # the cluster columns (e.g. already-wide count columns from prior pipeline)
        meta = {'feature', 'Feature', 'Trait', 'trait', 'Loci', 'loci',
                'Locus', 'locus', 'SNP', 'snp', 'snp1',
                'snp1_gene_names', 'snp2_gene_names',
                'snp1_gene_name',  'snp2_gene_name'}
        cluster_cols = [c for c in df.columns if c not in meta]
        if not cluster_cols:
            file_stem = os.path.splitext(os.path.basename(filepath))[0]
            df[file_stem] = 1
            cluster_cols = [file_stem]
        return df[[c for c in id_cols + cluster_cols if c in df.columns]], \
               cluster_cols, feat_col, loci_col, trait_col


# =============================================================================
#  STEP 5: BUILD CLUSTER LOOKUP AND MAP FEATURES
# =============================================================================

def _norm(s):
    """
    Thin wrapper: apply full canonical normalisation to a single value.
    Returns empty string for null/None.
    """
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ''
    return normalise_name(str(s))


def build_lookup(wide_df, cluster_cols, feat_col, loci_col, trait_col):
    """
    Build per-cluster lookup dicts of normalised identifier strings.

    All values are stored in their canonical form (via normalise_name) so
    that spelling variants collapse to the same key at query time.

    Returns
    -------
    dict  {cluster_col: {'features': set, 'loci': set}}
         'features' – values from the feature/SNP identifier column
         'loci'     – values from the Loci/loci column
    """
    lookup = {}
    for cc in cluster_cols:
        if cc not in wide_df.columns:
            continue
        # Rows where this cluster is represented (non-null, non-zero)
        present = wide_df[wide_df[cc].notna() & (wide_df[cc] != 0)]

        features_set, loci_set = set(), set()

        if feat_col and feat_col in present.columns:
            for v in present[feat_col].dropna():
                # Split on commas/semicolons only — NOT whitespace.
                # rsIDs never contain spaces; clinical names do ("Fasting glucose",
                # "Hemoglobin concentration", etc.) and must be kept whole so that
                # normalise_name() can convert the full phrase to a single canonical key
                # (e.g. "Hemoglobin concentration" → "haemoglobin_concentration").
                for token in re.split(r'[,;]+', str(v)):
                    token = token.strip()
                    n = _norm(token)
                    if n:
                        features_set.add(n)

        if loci_col and loci_col in present.columns:
            for v in present[loci_col].dropna():
                n = _norm(v)
                if n:
                    loci_set.add(n)

        lookup[cc] = {'features': features_set, 'loci': loci_set}

    return lookup


def map_row_to_clusters(feature, snp1, snp2, snp1_genes, snp2_genes,
                         feature_source, all_lookups):
    """
    Map a single feature to all functional clusters with a strict match order.

    Our side                         Cluster-file side
    ────────────────────────────     ─────────────────────────────────────────
    snp1  (primary SNP / index SNP)  feature col  (rsID or clinical measure,
    snp2  (secondary SNP, GxG)       ──────────────  renamed from Trait where needed)
    feature (raw name, clinical      loci col     (gene names; equivalent to our
             fallback only)          ──────────────  snp1_gene_names / snp2_gene_names)

    Match order — first hit wins per cluster column:
      1. feature → feature : snp1 (then snp2) vs cluster 'feature' col.
                             For clinical rows with no SNP annotation, raw 'feature'
                             name is used instead (direct trait-name match).
      2. loci    → loci    : snp1_gene_names + snp2_gene_names tokens (our loci)
                             vs cluster 'loci' col (gene-name loci on cluster side).
      3. feature → loci    : snp1 / snp2 rsID tokens vs cluster 'loci' col
                             (fallback: SNP may not be the index SNP but shares a locus).

    Parameters
    ----------
    all_lookups : list of (col_key, lookup_dict) where lookup_dict has keys
                  'features', 'loci'

    Returns
    -------
    hits        : dict  {col_key: 0 or 1}
    method_used : str   'feature_feature' | 'loci_loci' | 'feature_loci' | 'none'
                  Reflects the most specific method that produced any positive hit
                  across all cluster columns for this row.
    """
    # --- Our "feature" candidates (step 1 and step 3) -----------------------
    # Primary: snp1, secondary: snp2.  Only use raw feature name when no SNP
    # annotation exists (i.e. clinical measures).
    snp_candidates = set()
    for v in [snp1, snp2]:
        n = _norm(v)
        if n:
            snp_candidates.add(n)

    # Clinical fallback: feature name used for feature→feature when no SNPs
    clinical_candidates = set()
    if not snp_candidates:
        n = _norm(feature)
        if n:
            clinical_candidates.add(n)

    feature_candidates = snp_candidates | clinical_candidates

    # --- Our "loci" candidates (step 2) -------------------------------------
    # snp1_gene_names and snp2_gene_names are our loci — semicolon-separated
    # gene names that map to the 'loci' column on the cluster-file side.
    loci_candidates = set()
    for gene_str in [snp1_genes, snp2_genes]:
        if gene_str and pd.notna(gene_str):
            for g in re.split(r'[;,]', str(gene_str)):
                n = _norm(g.strip())
                if n:
                    loci_candidates.add(n)

    # --- Map each cluster column ---------------------------------------------
    hits = {}
    best_method = 'none'
    method_rank = {'none': 0, 'feature_loci': 1, 'loci_loci': 2, 'feature_feature': 3}

    for col_key, lookup in all_lookups:
        hit = 0
        method = 'none'

        # 1. feature → feature
        #    snp1/snp2 rsID (or clinical name) vs cluster feature col
        if feature_candidates & lookup['features']:
            hit = 1
            method = 'feature_feature'

        # 2. loci → loci
        #    our gene-name tokens (snp1_gene_names / snp2_gene_names)
        #    vs cluster loci col (gene-name loci on cluster-file side)
        elif loci_candidates & lookup['loci']:
            hit = 1
            method = 'loci_loci'

        # 3. feature → loci
        #    snp1/snp2 rsID vs cluster loci col
        #    (SNP not listed as index SNP but appears in a locus entry)
        elif snp_candidates & lookup['loci']:
            hit = 1
            method = 'feature_loci'

        hits[col_key] = hit
        if method_rank[method] > method_rank[best_method]:
            best_method = method

    return hits, best_method


# =============================================================================
#  MAIN
# =============================================================================

def main(pheno_data, raw_features_file,
         filter_strategy='tiered',
         use_set='holdout',
         balance_strategy='unbalanced',
         n_iterations=1000,
         alpha=0.05,
         min_effect_size=0.1,
         use_all=False):

    scores_path = os.path.join(pheno_data, 'scores')
    output_path = os.path.join(scores_path, 'functionalEnrichment')
    os.makedirs(output_path, exist_ok=True)

    # -------------------------------------------------------------------------
    # Load PRS bin assignments
    # -------------------------------------------------------------------------
    if use_set == 'holdout':
        prs_file = os.path.join(scores_path, 'combinedPRSGroups.holdout.filtered.csv')
        if not os.path.exists(prs_file):
            prs_file = os.path.join(scores_path, 'combinedPRSGroups.filtered.csv')
    else:
        prs_file = os.path.join(scores_path, 'combinedPRSGroups.filtered.csv')

    print(f"\nLoading PRS bins: {prs_file}")
    combined_prs = pd.read_csv(prs_file)

    bin_cols = [c for c in combined_prs.columns
                if c.startswith('bin_') and 'combined' not in c]
    if not use_all:
        bin_cols = [c for c in bin_cols if 'all' not in c]

    combined_prs = combined_prs[['IID', 'PHENOTYPE'] + bin_cols].set_index('IID')
    cohort_names = [c.replace('bin_', '') for c in bin_cols]
    print(f"Cohorts: {cohort_names}")

    # -------------------------------------------------------------------------
    # Identify exclusive high-risk cases + controls per cohort
    # -------------------------------------------------------------------------
    print("\n--- Step 1: Identifying exclusive high-risk cases ---")
    cohort_groups = identify_exclusive_groups(combined_prs, bin_cols)

    # -------------------------------------------------------------------------
    # Load raw clinical features
    # -------------------------------------------------------------------------
    print(f"\nLoading raw features: {raw_features_file}")
    raw_df = pd.read_csv(raw_features_file)
    if 'Participant ID' in raw_df.columns:
        raw_df.rename(columns={'Participant ID': 'IID'}, inplace=True)
    raw_df.set_index('IID', inplace=True)

    exclude_meta = {'PHENOTYPE', 'SEX', 'age'} | set(bin_cols) | set(cohort_names)
    feature_cols_raw = [c for c in raw_df.columns if c not in exclude_meta]

    # Normalise clinical feature column names to collapse spelling variants.
    # e.g. "Residual Glycaemia", "residual_glycaemic", "residual glycemia" → same key.
    canon_map, canon_reverse = normalise_series_names(feature_cols_raw)
    # Rename raw_df columns using canonical names; keep first occurrence on collision.
    canonical_cols = {}
    rename_dict = {}
    for orig in feature_cols_raw:
        can = canon_map[orig]
        if can not in canonical_cols:
            canonical_cols[can] = orig
            rename_dict[orig] = can
        else:
            # Duplicate canonical name — drop the later occurrence, warn
            print(f"  [NORM] Dropping duplicate canonical '{can}': keeping "
                  f"'{canonical_cols[can]}', dropping '{orig}'")
    raw_df = raw_df[[c for c in raw_df.columns
                     if c in exclude_meta or c in rename_dict]]
    raw_df.rename(columns=rename_dict, inplace=True)
    feature_cols = [c for c in raw_df.columns if c not in exclude_meta]
    print(f"Clinical feature columns: {len(feature_cols_raw)} raw → "
          f"{len(feature_cols)} after normalisation")

    # -------------------------------------------------------------------------
    # Load pre-computed genomic features
    # -------------------------------------------------------------------------
    print("\n--- Step 2b: Loading pre-computed genomic SHAP features ---")
    genomic_df = load_genomic_features(scores_path, filter_strategy)

    # -------------------------------------------------------------------------
    # Load feature annotation
    # -------------------------------------------------------------------------
    print("\n--- Step 3: Loading feature annotation ---")
    annotation_df = load_feature_annotation(scores_path)

    # -------------------------------------------------------------------------
    # Find important features per cohort
    # -------------------------------------------------------------------------
    print(f"\n--- Step 2a: Finding important clinical features "
          f"(strategy: {balance_strategy}) ---")

    all_feature_rows = []

    for cohort, groups in cohort_groups.items():
        cases = [i for i in groups['cases'] if i in raw_df.index]
        controls = [i for i in groups['controls'] if i in raw_df.index]

        print(f"\n  {cohort}: {len(cases)} cases | {len(controls)} controls "
              f"(after raw-feature filter)")

        # ---- Clinical features ----
        if len(cases) >= 5 and len(controls) >= 5:
            if balance_strategy == 'balanced':
                clin_result = find_clinical_features_balanced(
                    raw_df, cases, controls, feature_cols,
                    alpha=alpha, min_effect_size=min_effect_size,
                    n_iterations=n_iterations,
                )
            else:
                clin_result = find_clinical_features_unbalanced(
                    raw_df, cases, controls, feature_cols,
                    alpha=alpha, min_effect_size=min_effect_size,
                )

            if not clin_result.empty:
                sig_clin = clin_result[clin_result['significant']].copy()
                print(f"    Clinical: {len(sig_clin)}/{len(clin_result)} significant")
                sig_clin['cohort'] = cohort
                sig_clin['raw_feature'] = sig_clin['feature']
                sig_clin['feature_source'] = 'clinical'
                sig_clin['odds_ratio'] = np.nan
                sig_clin['specificity_tier'] = np.nan
                all_feature_rows.append(sig_clin)
        else:
            print(f"    [SKIP clinical] too few samples")

        # ---- Genomic features (pre-computed SHAP, cohort-specific) ----
        if not genomic_df.empty:
            cohort_genomic = genomic_df[genomic_df['cohort'] == cohort].copy()
            if not cohort_genomic.empty:
                cohort_genomic['raw_feature'] = cohort_genomic['feature']
                cohort_genomic['feature_source'] = 'genomic_shap'
                cohort_genomic['p_value'] = np.nan
                cohort_genomic['p_value_fdr'] = np.nan
                cohort_genomic['effect_size_r'] = np.nan
                cohort_genomic['significant'] = True
                print(f"    Genomic SHAP: {len(cohort_genomic)} features")
                all_feature_rows.append(cohort_genomic)

    if not all_feature_rows:
        print("\n[ERROR] No important features found. Check inputs and thresholds.")
        return

    # -------------------------------------------------------------------------
    # Combine, deduplicate, and expand to one row per SNP
    # -------------------------------------------------------------------------
    keep_cols = ['raw_feature', 'feature', 'cohort', 'feature_source',
                 'p_value', 'p_value_fdr', 'effect_size_r',
                 'odds_ratio', 'specificity_tier']

    combined_features = []
    for df in all_feature_rows:
        present = [c for c in keep_cols if c in df.columns]
        sub = df[present].copy()
        for c in keep_cols:
            if c not in sub.columns:
                sub[c] = np.nan
        combined_features.append(sub[keep_cols])

    feature_df = pd.concat(combined_features, ignore_index=True)
    # Deduplicate same raw_feature × cohort × source
    feature_df.drop_duplicates(
        subset=['raw_feature', 'cohort', 'feature_source'], inplace=True
    )
    print(f"\nTotal features (before SNP expansion): {len(feature_df)}")

    # -------------------------------------------------------------------------
    # Expand GxG features and merge gene annotation
    # -------------------------------------------------------------------------
    print("--- Step 3: Expanding to individual SNPs and merging annotation ---")
    expanded_df = expand_to_snp_rows(feature_df, annotation_df)
    # Deduplicate expanded (same snp can appear in snp1 of one row and snp2 of another)
    expanded_df.drop_duplicates(
        subset=['feature', 'cohort', 'feature_source'], inplace=True
    )
    print(f"Features after SNP expansion: {len(expanded_df)}")

    # -------------------------------------------------------------------------
    # Load functional cluster files
    # -------------------------------------------------------------------------
    print("\n--- Step 4: Loading functional cluster files ---")
    subtype_path = os.path.join(os.path.dirname(pheno_data), 'subtypeResearch')

    # Populate clinical alias map BEFORE any normalise_name calls on cluster data
    global _CLINICAL_FEATURE_MAP
    _CLINICAL_FEATURE_MAP = load_clinical_feature_map(subtype_path)

    cluster_pattern = os.path.join(subtype_path, '*functional_cluster*.csv')
    cluster_files = sorted(glob.glob(cluster_pattern))

    if not cluster_files:
        print(f"  [WARN] No cluster files found: {cluster_pattern}")

    # flat_lookups: list of (col_key, lookup_dict) for map_row_to_clusters
    flat_lookups = []
    all_cluster_col_keys = []

    for fpath in cluster_files:
        fname = os.path.basename(fpath)
        file_stem = re.sub(r'[^a-zA-Z0-9_]', '_', os.path.splitext(fname)[0])
        print(f"  {fname}")
        try:
            wide_df, cluster_cols, feat_col, loci_col, trait_col = load_cluster_file(fpath)
            lookup = build_lookup(wide_df, cluster_cols, feat_col, loci_col, trait_col)

            for cc in cluster_cols:
                # Prefix with file stem when multiple cluster files exist
                col_key = f"{file_stem}__{cc}" if len(cluster_files) > 1 else cc
                if col_key not in all_cluster_col_keys:
                    all_cluster_col_keys.append(col_key)
                if cc in lookup:
                    flat_lookups.append((col_key, lookup[cc]))

            n_features = sum(len(v['features']) for v in lookup.values())
            print(f"    → {len(cluster_cols)} clusters, {n_features} feature entries, "
                  f"{sum(len(v['loci']) for v in lookup.values())} loci entries")
        except Exception as exc:
            print(f"  [ERROR] {fname}: {exc}")

    # -------------------------------------------------------------------------
    # Map features to clusters
    # -------------------------------------------------------------------------
    print("\n--- Step 5: Mapping features to functional clusters ---")
    output_rows = []
    for _, row in expanded_df.iterrows():
        feat_val = row.get('feature', row.get('raw_feature'))
        cluster_hits, match_method = map_row_to_clusters(
            feature=feat_val,
            snp1=row.get('snp1'),
            snp2=row.get('snp2'),
            snp1_genes=row.get('snp1_gene_names'),
            snp2_genes=row.get('snp2_gene_names'),
            feature_source=row.get('feature_source'),
            all_lookups=flat_lookups,
        )
        # Count how many functional cluster columns have a positive hit (1)
        n_matched = sum(cluster_hits.values())

        output_rows.append({
            'raw_feature':        row.get('raw_feature', feat_val),
            'feature':            feat_val,
            'cohort':             row['cohort'],
            'feature_source':     row.get('feature_source'),
            'p_value':            row.get('p_value'),
            'p_value_fdr':        row.get('p_value_fdr'),
            'effect_size_r':      row.get('effect_size_r'),
            'odds_ratio':         row.get('odds_ratio'),
            'specificity_tier':   row.get('specificity_tier'),
            'snp1_gene_names':    row.get('snp1_gene_names'),
            'snp2_gene_names':    row.get('snp2_gene_names'),
            'match_method':       match_method,
            'n_clusters_matched': n_matched,
            **cluster_hits,        # binary 0/1 per functional cluster column
        })

    if not output_rows:
        print("[ERROR] No output rows generated.")
        return

    output_df = pd.DataFrame(output_rows)

    # Ensure every cluster column is present (fill 0 if a file errored out)
    for ck in all_cluster_col_keys:
        if ck not in output_df.columns:
            output_df[ck] = 0

    # Enforce binary integers (guard against float bleed from pivot_table).
    # A non-null, non-zero value → 1; null or zero → 0.
    for ck in all_cluster_col_keys:
        output_df[ck] = (output_df[ck].fillna(0) != 0).astype(int)

    # Recompute n_clusters_matched after enforcing binary ints
    if all_cluster_col_keys:
        output_df['n_clusters_matched'] = output_df[all_cluster_col_keys].sum(axis=1)

    # Canonical column order: metadata | match tracking | cluster binary columns
    meta_cols = [
        'raw_feature', 'feature', 'cohort', 'feature_source',
        'p_value', 'p_value_fdr', 'effect_size_r', 'odds_ratio',
        'specificity_tier', 'snp1_gene_names', 'snp2_gene_names',
        'match_method', 'n_clusters_matched',
    ]
    cluster_cols_out = [c for c in output_df.columns if c not in meta_cols]
    output_df = output_df[meta_cols + cluster_cols_out]

    # -------------------------------------------------------------------------
    # Save and summarise
    # -------------------------------------------------------------------------
    out_file = os.path.join(
        output_path,
        f'functional_cluster_enrichment_{use_set}_{balance_strategy}.csv'
    )
    output_df.to_csv(out_file, index=False)

    n_any_cluster = int((output_df['n_clusters_matched'] > 0).sum()) \
        if 'n_clusters_matched' in output_df.columns else 0
    print(f"\nSaved → {out_file}")
    print(f"Shape : {output_df.shape}")
    print(f"Features mapped to ≥1 cluster: {n_any_cluster}/{len(output_df)}")

    print("\nPer-cohort feature counts (cohort × source):")
    summary = (
        output_df.groupby(['cohort', 'feature_source'])
        .size()
        .rename('n_features')
        .reset_index()
    )
    print(summary.to_string(index=False))

    print("\nMatch method breakdown:")
    method_summary = (
        output_df.groupby('match_method')
        .size()
        .rename('n_features')
        .reset_index()
    )
    print(method_summary.to_string(index=False))

    if cluster_cols_out:
        print("\nPer-cluster mapping counts (features with binary=1):")
        counts = output_df[cluster_cols_out].sum().astype(int).reset_index()
        counts.columns = ['cluster', 'n_features_mapped']
        counts = counts[counts['n_features_mapped'] > 0].sort_values(
            'n_features_mapped', ascending=False
        )
        if not counts.empty:
            print(counts.to_string(index=False))
        else:
            print("  (no features mapped to any cluster)")

    return output_df


# =============================================================================
#  CLI
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Enrich cohort-exclusive high-risk case features against "
            "functional cluster annotations."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--pheno_data',
        help='Path to combinedAnalysis directory (e.g. results/type2Diabetes/combinedAnalysis); '
             'subtypeResearch/ is expected one level up (results/type2Diabetes/subtypeResearch/)',
    )
    parser.add_argument(
        '--raw_features_file',
        help='Path to raw (untransformed) clinical features CSV',
    )
    parser.add_argument(
        '--filter_strategy', default='tiered',
        choices=['strict', 'majority', 'differential', 'tiered'],
        help='Strategy used when loading pre-computed SHAP genomic features (default: tiered)',
    )
    parser.add_argument(
        '--use_set', default='holdout',
        choices=['holdout', 'validation', 'all'],
        help='Which individual set to use for clinical comparisons (default: holdout)',
    )
    parser.add_argument(
        '--balance_strategy', default='unbalanced',
        choices=['unbalanced', 'balanced'],
        help=(
            'unbalanced (default): all controls vs exclusive cases — '
            'maximises power when control pool >> n_cases. '
            'balanced: bootstrap n_cases controls over n_iterations — '
            'recommended when pool is only 2–5× cases.'
        ),
    )
    parser.add_argument(
        '--n_iterations', type=int, default=1000,
        help='Number of bootstrap iterations for balanced strategy (default: 1000)',
    )
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='FDR significance threshold (default: 0.05)',
    )
    parser.add_argument(
        '--min_effect_size', type=float, default=0.1,
        help='Minimum |rank-biserial r| to flag as meaningful (default: 0.1)',
    )
    parser.add_argument(
        '--use_all', action='store_true',
        help='Include "all" cohort bins (excluded by default)',
    )

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get('PHENO_DATA')
    raw_features_file = args.raw_features_file or os.environ.get('RAW_FEATURES_FILE')

    if not pheno_data:
        parser.error('Provide --pheno_data or set PHENO_DATA env var')
    if not raw_features_file:
        parser.error('Provide --raw_features_file or set RAW_FEATURES_FILE env var')

    print(f"pheno_data       : {pheno_data}")
    print(f"raw_features_file: {raw_features_file}")
    print(f"filter_strategy  : {args.filter_strategy}")
    print(f"use_set          : {args.use_set}")
    print(f"balance_strategy : {args.balance_strategy}")
    print(f"n_iterations     : {args.n_iterations}")
    print(f"alpha            : {args.alpha}")
    print(f"min_effect_size  : {args.min_effect_size}")
    print(f"use_all          : {args.use_all}")

    main(
        pheno_data=pheno_data,
        raw_features_file=raw_features_file,
        filter_strategy=args.filter_strategy,
        use_set=args.use_set,
        balance_strategy=args.balance_strategy,
        n_iterations=args.n_iterations,
        alpha=args.alpha,
        min_effect_size=args.min_effect_size,
        use_all=args.use_all,
    )
