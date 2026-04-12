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
       Clinical features : Mann-Whitney U on raw_feature_scores
                           (unbalanced = all controls; balanced = bootstrap n_cases
                            controls over --n_iterations iterations)
       Genomic features  : loaded from pre-computed SHAP output
                           (importantFeaturesAcrossCohortsAndTrainingData.Filtered.*.csv)
                           — already cohort-exclusive by cross-cohort filtering strategy

  3. Merge genomic features with feature_scores_annotated.csv to get snp1/snp2
     and their gene names. GxG features (two rsIDs) expand to one row per SNP.

  4. Load functional cluster CSVs from {results_path}/subtypeResearch/:
       - all *.csv files except clinical_feature_map.csv and mahajan_clinical_loci_map.csv
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
    --raw_feature_scores /path/to/raw_features.csv \\
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

_RSID_RE = re.compile(r'^rs\d+$', re.IGNORECASE)


def _infer_feature_source(feature: str) -> str:
    """Return 'genomic_shap' or 'clinical' from a feature name.

    Rules (in order):
      • All comma-separated tokens are rsIDs  → genomic_shap  (G or GxG)
      • Single rsID                           → genomic_shap
      • Two-part underscore token (HLA allele)→ genomic_shap
      • Anything else                         → clinical
    """
    feature = str(feature).strip()
    parts = [p.strip() for p in feature.split(',')]
    if all(_RSID_RE.match(p) for p in parts if p):
        return 'genomic_shap'
    if len(parts) == 1:
        # HLA allele: two alphanumeric tokens joined by '_', neither is an rsID
        hla_parts = feature.split('_')
        if (len(hla_parts) == 2
                and all(p and p.isalnum() for p in hla_parts)
                and not _RSID_RE.match(feature)):
            return 'genomic_shap'
    return 'clinical'


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


def filter_clinical_features_cross_cohort(
    cohort_results: dict,
    majority_threshold: float = 0.6,
    delta_threshold: float = 0.1,
) -> pd.DataFrame:
    """
    Apply cross-cohort specificity filtering to per-cohort MWU results.

    Mirrors the tiered PRS strategy but uses effect_size_r instead of z_score:
      Tier 1 (strict)      – significant in matching cohort, NOT significant in
                             any other cohort.
      Tier 2 (majority)    – significant in matching cohort, NOT significant in
                             >= majority_threshold fraction of other cohorts.
      Tier 3 (differential)– matching cohort has the highest |effect_size_r| AND
                             that value minus the mean |r| of other cohorts >=
                             delta_threshold.

    For each feature the tier is determined and it is assigned to the single
    cohort with the highest |effect_size_r| (among tiers 1/2/3 candidates).

    Parameters
    ----------
    cohort_results : dict  {cohort_name: result_df}
        Each result_df has columns: feature, p_value_fdr, effect_size_r,
        significant (bool).
    majority_threshold : float  fraction of other cohorts that must be
        non-significant for tier 2.
    delta_threshold : float  minimum |r| advantage over mean-of-others for tier 3.

    Returns
    -------
    DataFrame with columns: feature, cohort, effect_size_r, p_value_fdr,
                             specificity_tier  (1 / 2 / 3),
                             n_cohorts, feature_source='clinical',
                             training_data_used='Clinical'
    """
    cohorts = list(cohort_results.keys())
    if not cohorts:
        return pd.DataFrame()

    # Collect all features that were significant in at least one cohort
    all_features = set()
    for df in cohort_results.values():
        if df.empty:
            continue
        all_features.update(df.loc[df['significant'], 'feature'].tolist())

    if not all_features:
        print("  [clinical cross-cohort] No significant features in any cohort.")
        return pd.DataFrame()

    # Build lookup: (cohort, feature) → row dict for fast access
    result_lookup = {}
    for cohort, df in cohort_results.items():
        if df.empty:
            continue
        for _, row in df.iterrows():
            result_lookup[(cohort, row['feature'])] = row

    output_rows = []
    for feat in all_features:
        # Collect per-cohort effect sizes and significance flags
        cohort_r    = {}   # cohort → |effect_size_r|  (NaN if not tested)
        cohort_sig  = {}   # cohort → bool significant
        cohort_fdr  = {}   # cohort → p_value_fdr

        for cohort in cohorts:
            row = result_lookup.get((cohort, feat))
            if row is not None:
                cohort_r[cohort]   = abs(float(row['effect_size_r'])) if pd.notna(row.get('effect_size_r')) else np.nan
                cohort_sig[cohort] = bool(row.get('significant', False))
                cohort_fdr[cohort] = float(row['p_value_fdr']) if pd.notna(row.get('p_value_fdr')) else 1.0
            else:
                cohort_r[cohort]   = np.nan
                cohort_sig[cohort] = False
                cohort_fdr[cohort] = 1.0

        sig_cohorts  = [c for c in cohorts if cohort_sig[c]]
        n_sig        = len(sig_cohorts)
        n_cohorts    = len(cohorts)

        if n_sig == 0:
            continue

        # Determine best cohort: highest |r| among significant cohorts
        best_cohort = max(sig_cohorts, key=lambda c: cohort_r[c] if not np.isnan(cohort_r[c]) else -1)
        best_r      = cohort_r[best_cohort]
        other_r     = [cohort_r[c] for c in cohorts if c != best_cohort and not np.isnan(cohort_r[c])]
        n_other_sig = sum(1 for c in cohorts if c != best_cohort and cohort_sig[c])

        tier = None

        # Tier 1: significant ONLY in best_cohort
        if n_other_sig == 0:
            tier = 1

        # Tier 2: not significant in >= majority_threshold of other cohorts
        if tier is None and n_cohorts > 1:
            frac_other_non_sig = (n_cohorts - 1 - n_other_sig) / (n_cohorts - 1)
            if frac_other_non_sig >= majority_threshold:
                tier = 2

        # Tier 3: highest |r| AND delta above mean of others >= delta_threshold
        if tier is None and other_r:
            mean_other_r = float(np.nanmean(other_r))
            if not np.isnan(best_r) and (best_r - mean_other_r) >= delta_threshold:
                tier = 3

        if tier is None:
            continue  # feature not specific enough to any cohort

        best_row = result_lookup.get((best_cohort, feat), {})
        output_rows.append({
            'feature':            feat,
            'raw_feature':        feat,
            'cohort':             best_cohort,
            'effect_size_r':      float(best_row.get('effect_size_r', np.nan)),
            'p_value':            float(best_row.get('p_value', np.nan)),
            'p_value_fdr':        cohort_fdr[best_cohort],
            'specificity_tier':   tier,
            'n_cohorts':          n_cohorts,
            'feature_source':     'clinical',
            'training_data_used': 'Clinical',
            'odds_ratio':         np.nan,
            'n_exclusive_cohorts': np.nan,
            'exclusive_cohorts':   np.nan,
        })

    if not output_rows:
        return pd.DataFrame()

    result = pd.DataFrame(output_rows)
    print(f"  Clinical cross-cohort filter: {len(all_features)} features tested "
          f"→ {len(result)} survived")
    tier_counts = result['specificity_tier'].value_counts().sort_index()
    for t, n in tier_counts.items():
        print(f"    Tier {t}: {n} features")
    return result


# =============================================================================
#  STEP 2b: GENOMIC FEATURES (from pre-computed SHAP)
# =============================================================================

def load_genomic_features(scores_path, filter_strategy):
    """
    Load SHAP-based cohort-exclusive genomic features from pre-computed output.
    These are already filtered for cross-cohort exclusivity.

    Returns DataFrame with columns: feature, cohort, feature_source,
                                     odds_ratio, [specificity_tier],
                                     [training_data_used],
                                     [n_exclusive_cohorts], [exclusive_cohorts]
    """
    path = os.path.join(
        scores_path, 'importantCohortScores',
        f'importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv'
    )
    if not os.path.exists(path):
        print(f"  [WARN] Filtered features file not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    # Normalise cohort column.
    # HighCases/LowControls rows come from find_tiered_features which writes
    # 'feature_set' (= the matching cohort name) but no 'cohort' column.
    # ExclusiveHighCases rows come from run_exclusive_high_bootstrap which writes
    # 'cohort' but no 'feature_set'.
    # When both training categories are concatenated in the CSV, both columns
    # exist but each set of rows has NaN in the column belonging to the other.
    # Fill 'cohort' from 'feature_set' wherever cohort is missing, then drop
    # the now-redundant 'feature_set' column.
    if 'feature_set' in df.columns:
        if 'cohort' not in df.columns:
            df.rename(columns={'feature_set': 'cohort'}, inplace=True)
        else:
            df['cohort'] = df['cohort'].fillna(df['feature_set'])
            df.drop(columns=['feature_set'], inplace=True)

    # Include HighCases, LowControls, and ExclusiveHighCases
    if 'training_data_used' in df.columns:
        df = df[df['training_data_used'].isin(
            ['HighCases', 'LowControls', 'ExclusiveHighCases']
        )]

    # Drop any rows still missing a cohort (shouldn't happen, but guard)
    n_before = len(df)
    df = df[df['cohort'].notna()]
    if len(df) < n_before:
        print(f"  [WARN] Dropped {n_before - len(df)} rows with missing cohort")

    # ── Issue 3: ExclusiveHighCases odds_ratio ────────────────────────────────
    # ExclusiveHighCases rows are produced by a different bootstrap procedure
    # that writes coefs (log-odds coefficients) but does not pre-compute
    # odds_ratio.  Derive it from coefs so downstream ranking works correctly.
    if 'coefs' in df.columns and 'odds_ratio' in df.columns:
        exc_mask = (
            df.get('training_data_used', pd.Series()) == 'ExclusiveHighCases'
        ) & df['odds_ratio'].isna() & df['coefs'].notna()
        if exc_mask.any():
            df.loc[exc_mask, 'odds_ratio'] = np.exp(df.loc[exc_mask, 'coefs'])
            print(f"  Computed odds_ratio from coefs for "
                  f"{exc_mask.sum()} ExclusiveHighCases rows")

    # ── Issue 1: LowControls direction ───────────────────────────────────────
    # LowControls features represent the protective signature of each cohort
    # (bottom 20 % PRS).  Only features with OR < 1 are biologically meaningful
    # as protective associations; OR > 1 features from this stratum are risk-
    # allele carry-overs from the shared classifier and should not be reported
    # as protective evidence.
    if 'training_data_used' in df.columns and 'odds_ratio' in df.columns:
        lc_risk_mask = (
            (df['training_data_used'] == 'LowControls')
            & df['odds_ratio'].notna()
            & (df['odds_ratio'] >= 1.0)
        )
        n_lc_risk = lc_risk_mask.sum()
        if n_lc_risk:
            df = df[~lc_risk_mask].copy()
            print(f"  Removed {n_lc_risk} LowControls features with OR ≥ 1 "
                  f"(risk direction, not protective)")

    # Infer feature_source from the feature name itself
    df['feature_source'] = df['feature'].apply(_infer_feature_source)

    print(f"  Loaded {len(df)} filtered features ({df['cohort'].nunique()} cohorts)")
    if 'training_data_used' in df.columns:
        print('  ' + df.groupby('training_data_used').size()
              .rename('n').reset_index().to_string(index=False))
    return df


# =============================================================================
#  STEP 3: MERGE WITH feature_scores_annotated.csv AND EXPAND SNPs
# =============================================================================

def load_feature_annotation(scores_path):
    """Load feature_scores_annotated.csv.

    Searches in order:
      1. {scores_path}/feature_scores_annotated.csv  (e.g. combinedAnalysis/scores/)
      2. {pheno_data}/feature_scores_annotated.csv   (parent = combinedAnalysis/)
      3. {type2Diabetes}/feature_scores_annotated.csv (grandparent)
    """
    candidates = [
        os.path.join(scores_path, 'feature_scores_annotated.csv'),
        os.path.join(os.path.dirname(scores_path), 'feature_scores_annotated.csv'),
        os.path.join(os.path.dirname(os.path.dirname(scores_path)),
                     'feature_scores_annotated.csv'),
    ]
    for path in candidates:
        if os.path.exists(path):
            df = pd.read_csv(path)
            print(f"  Loaded feature_scores_annotated: {len(df)} rows  [{path}]")
            return df
    print(f"  [WARN] feature_scores_annotated.csv not found. Searched:\n"
          + '\n'.join(f"    {p}" for p in candidates))
    return pd.DataFrame()


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
    # Build fast lookup: normalised feature string → annotation row
    # Keys are strip()+lower() so minor case/whitespace differences don't cause misses.
    ann_lookup = {}
    if not annotation_df.empty and 'feature' in annotation_df.columns:
        for _, arow in annotation_df.iterrows():
            ann_lookup[str(arow['feature']).strip().lower()] = arow

    expanded = []
    n_annotated = 0
    n_unannotated_genomic = 0
    for _, row in feature_df.iterrows():
        raw = str(row['raw_feature'])
        # raw_feature may be NaN when load_genomic_features provided 'feature' only;
        # in that case fall back to the 'feature' column for annotation lookup.
        if raw.strip().lower() == 'nan' or not raw.strip():
            raw = str(row.get('feature', ''))
        base = row.to_dict()
        ann = ann_lookup.get(raw.strip().lower())

        snp1 = ann.get('snp1') if ann is not None and pd.notna(ann.get('snp1')) else None
        snp2 = ann.get('snp2') if ann is not None and pd.notna(ann.get('snp2')) else None
        snp1_genes = ann.get('snp1_gene_names') if ann is not None else None
        snp2_genes = ann.get('snp2_gene_names') if ann is not None else None

        if ann is not None:
            n_annotated += 1
        elif row.get('feature_source') == 'genomic_shap':
            n_unannotated_genomic += 1

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

    print(f"  Annotation lookup: {n_annotated}/{len(feature_df)} features matched in "
          f"feature_scores_annotated.csv")
    if n_unannotated_genomic:
        print(f"  [WARN] {n_unannotated_genomic} genomic features had no annotation entry "
              f"— loci-based cluster matching will be unavailable for those.")

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

    Assignment rule — max-cluster wins
    ------------------------------------
    Each row (feature / locus) is assigned to the ONE cluster column whose
    absolute beta value is the largest.  This prevents a trait like
    "waist circumference" (non-zero in every cluster but dominant in obesity)
    from matching every cluster equally.

    When a file has only a single cluster column (or the column is binary
    0/1 rather than continuous betas), every non-zero row is assigned to
    that column as before.

    All values are stored in their canonical form (via normalise_name) so
    that spelling variants collapse to the same key at query time.

    Returns
    -------
    dict  {cluster_col: {'features': set, 'loci': set}}
         'features' – values from the feature/SNP identifier column
         'loci'     – values from the Loci/loci column
    """
    # Initialise empty sets for every cluster column
    lookup = {cc: {'features': set(), 'loci': set()}
              for cc in cluster_cols if cc in wide_df.columns}

    # Only consider cluster columns that actually exist in the dataframe
    valid_cols = [cc for cc in cluster_cols if cc in wide_df.columns]
    if not valid_cols:
        return lookup

    # -------------------------------------------------------------------------
    # Determine the dominant cluster for each row.
    # For rows where all cluster values are zero / NaN → skip (no assignment).
    # -------------------------------------------------------------------------
    beta_df = wide_df[valid_cols].copy()
    # Coerce to numeric; non-numeric (e.g. text) → NaN
    beta_df = beta_df.apply(pd.to_numeric, errors='coerce')

    if len(valid_cols) == 1:
        # Single cluster: assign any non-zero/non-null row to it
        dominant = beta_df[valid_cols[0]].apply(
            lambda v: valid_cols[0] if pd.notna(v) and v != 0 else None
        )
    else:
        # Multi-cluster: argmax of absolute value per row
        abs_df = beta_df.abs()
        # idxmax returns NaN for all-NaN rows; replace by checking row max > 0
        dominant = abs_df.idxmax(axis=1)
        row_max = abs_df.max(axis=1)
        dominant = dominant.where(row_max > 0, other=None)

    # -------------------------------------------------------------------------
    # Populate lookup sets based on dominant cluster assignment
    # -------------------------------------------------------------------------
    for idx, dom_col in dominant.items():
        if dom_col is None:
            continue
        row = wide_df.loc[idx]

        if feat_col and feat_col in wide_df.columns:
            raw_feat = row.get(feat_col, None)
            if raw_feat is not None and not (isinstance(raw_feat, float) and np.isnan(raw_feat)):
                # Split on commas/semicolons only — NOT whitespace.
                # rsIDs never contain spaces; clinical names do ("Fasting glucose",
                # "Hemoglobin concentration") and must be kept whole so that
                # normalise_name() maps the full phrase to its canonical key.
                for token in re.split(r'[,;]+', str(raw_feat)):
                    n = _norm(token.strip())
                    if n:
                        lookup[dom_col]['features'].add(n)

        if loci_col and loci_col in wide_df.columns:
            raw_loci = row.get(loci_col, None)
            if raw_loci is not None and not (isinstance(raw_loci, float) and np.isnan(raw_loci)):
                n = _norm(str(raw_loci))
                if n:
                    lookup[dom_col]['loci'].add(n)

    return lookup


def load_mahajan_map(subtype_path: str) -> dict:
    """
    Load ``mahajan_clinical_loci_map.csv`` and build a lookup that maps each
    SNP rsID or gene-locus token to the clinical feature it is most strongly
    associated with (highest absolute beta across all clinical trait columns).

    CSV structure
    -------------
    feature  – rsID (e.g. rs2296172)
    loci     – gene symbol(s) (e.g. MACF1)
    <trait>  – one column per clinical trait; cells are numeric betas / effect
               sizes (NA allowed).

    Returns
    -------
    dict  {normalised_id: {'clinical_feature': normalised_trait_name,
                           'source': 'feature' | 'loci'}}

    ``source`` records whether the key came from the rsID (feature) column or
    the loci (gene) column so the downstream annotation can note the origin.
    """
    map_path = os.path.join(subtype_path, 'mahajan_clinical_loci_map.csv')
    if not os.path.isfile(map_path):
        print(f"  [WARN] mahajan_clinical_loci_map.csv not found at {map_path}; "
              "Mahajan clinical linkage disabled.")
        return {}

    df = pd.read_csv(map_path, encoding='utf-8-sig')
    # Drop completely empty columns (trailing Unnamed: cols)
    df = df.dropna(axis=1, how='all')
    df.columns = df.columns.str.strip()

    id_cols = {'feature', 'loci'}
    clinical_cols = [c for c in df.columns
                     if c not in id_cols and not c.startswith('Unnamed')]
    if not clinical_cols:
        print("  [WARN] No clinical columns found in mahajan_clinical_loci_map.csv.")
        return {}

    # Coerce clinical columns to numeric
    for c in clinical_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    mahajan: dict = {}
    n_mapped = 0

    for _, row in df.iterrows():
        # Find the clinical column with the maximum |beta| for this SNP/locus
        betas = row[clinical_cols].abs()
        if betas.isna().all():
            continue
        best_trait_col = betas.idxmax()
        best_trait_norm = normalise_name(best_trait_col)
        if not best_trait_norm:
            continue

        # Register the rsID key
        rsid = str(row.get('feature', '') or '').strip()
        if rsid and rsid.lower() not in ('nan', ''):
            key = _norm(rsid)
            if key:
                mahajan[key] = {'clinical_feature': best_trait_norm,
                                'source': 'feature'}
                n_mapped += 1

        # Register each gene token from the loci column
        loci_raw = str(row.get('loci', '') or '').strip()
        if loci_raw and loci_raw.lower() not in ('nan', ''):
            for gene in re.split(r'[;,\s]+', loci_raw):
                key = _norm(gene.strip())
                if key and key not in mahajan:
                    # Only add if not already registered by rsID (rsID takes priority)
                    mahajan[key] = {'clinical_feature': best_trait_norm,
                                    'source': 'loci'}
                    n_mapped += 1

    print(f"  Loaded mahajan_clinical_loci_map.csv: {len(df)} SNPs/loci, "
          f"{n_mapped} id→clinical mappings registered.")
    return mahajan


def apply_mahajan_clinical_features(expanded_df: pd.DataFrame,
                                     mahajan_map: dict) -> pd.DataFrame:
    """
    For each genomic row in *expanded_df*, look up snp1, snp2, and
    snp1_gene_names / snp2_gene_names tokens against *mahajan_map* and
    record the linked clinical feature.

    New columns added
    -----------------
    mahajan_clinical_feature : str   – canonical clinical feature name from
                                       Mahajan map (empty string if no hit)
    mahajan_source           : str   – 'snp1' | 'snp2' | 'gene' | '' –
                                       which identifier produced the hit

    Only genomic rows (feature_source == 'genomic_shap') are queried; clinical
    rows already have a clinical identity and are left unchanged.
    """
    mahajan_cf  = []
    mahajan_src = []

    for _, row in expanded_df.iterrows():
        cf, src = '', ''
        if row.get('feature_source') == 'genomic_shap' and mahajan_map:
            # Priority: snp1 → snp2 → gene tokens from snp1_genes → snp2_genes
            candidates = []
            for label, val in [('snp1', row.get('snp1')),
                                ('snp2', row.get('snp2'))]:
                n = _norm(val)
                if n:
                    candidates.append((label, n))
            for gene_col, label in [('snp1_gene_names', 'gene'),
                                     ('snp2_gene_names', 'gene')]:
                gene_str = row.get(gene_col, '')
                if gene_str and pd.notna(gene_str):
                    for g in re.split(r'[;,\s]+', str(gene_str)):
                        n = _norm(g.strip())
                        if n:
                            candidates.append((label, n))

            for label, key in candidates:
                if key in mahajan_map:
                    cf  = mahajan_map[key]['clinical_feature']
                    src = label
                    break

        mahajan_cf.append(cf)
        mahajan_src.append(src)

    expanded_df = expanded_df.copy()
    expanded_df['mahajan_clinical_feature'] = mahajan_cf
    expanded_df['mahajan_source']           = mahajan_src
    return expanded_df


def map_row_to_clusters(feature, snp1, snp2, snp1_genes, snp2_genes,
                         feature_source, all_lookups,
                         mahajan_clinical_feature=''):
    """
    Map a single feature to all functional clusters with a strict match order.

    Our side                         Cluster-file side
    ────────────────────────────     ─────────────────────────────────────────
    snp1  (primary SNP / index SNP)  feature col  (rsID or clinical measure,
    snp2  (secondary SNP, GxG)       ──────────────  renamed from Trait where needed)
    feature (raw name, clinical      loci col     (gene names; equivalent to our
             fallback only)          ──────────────  snp1_gene_names / snp2_gene_names)

    Match order — first hit wins per cluster column (most → least specific):
      1. feature → feature         : snp1/snp2 rsID (or clinical name for non-genomic
                                     rows) vs cluster 'feature' col.
      1b. feature → feature        : Mahajan-linked clinical trait for this SNP vs
           (mahajan)                 cluster 'feature' col.
      2. loci    → loci            : snp1_gene_names / snp2_gene_names gene tokens
                                     vs cluster 'loci' col.
      3. loci    → feature         : same gene tokens vs cluster 'feature' col
                                     (some cluster files use gene symbols as feature
                                     identifiers rather than rsIDs).
      4. feature → loci            : snp1/snp2 rsID tokens vs cluster 'loci' col
                                     (SNP not the index SNP but shares a locus).

    Parameters
    ----------
    all_lookups              : list of (col_key, lookup_dict) where lookup_dict has
                               keys 'features', 'loci'
    mahajan_clinical_feature : str   normalised clinical feature name derived from
                               mahajan_clinical_loci_map.csv for this SNP/locus;
                               empty string if none was found.

    Returns
    -------
    hits        : dict  {col_key: 0 or 1}
    method_used : str   'feature_feature' | 'feature_feature_mahajan' |
                        'loci_loci' | 'loci_feature' | 'feature_loci' | 'none'
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

    # Mahajan-derived clinical feature for this SNP/locus (genomic rows only).
    # Kept separate so we can distinguish the match method.
    mahajan_candidates = set()
    if mahajan_clinical_feature:
        n = _norm(mahajan_clinical_feature)
        if n:
            mahajan_candidates.add(n)

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
    method_rank = {
        'none':                     0,
        'feature_loci':             1,   # SNP rsID vs cluster loci col
        'loci_feature':             2,   # our gene names vs cluster feature col
        'loci_loci':                3,   # our gene names vs cluster loci col
        'feature_feature_mahajan':  4,   # clinical bridge via Mahajan map
        'feature_feature':          5,   # direct SNP / clinical name match (highest)
    }

    for col_key, lookup in all_lookups:
        hit = 0
        method = 'none'

        # 1. feature → feature (direct)
        #    snp1/snp2 rsID (or clinical name) vs cluster feature col
        if feature_candidates & lookup['features']:
            hit = 1
            method = 'feature_feature'

        # 1b. feature → feature via Mahajan clinical bridge
        #     SNP linked to a clinical trait in Mahajan et al. → that trait
        #     is matched against the cluster feature col.  Ranked below direct
        #     SNP / name match but above loci-level matching.
        elif mahajan_candidates & lookup['features']:
            hit = 1
            method = 'feature_feature_mahajan'

        # 2. loci → loci
        #    our gene-name tokens (snp1_gene_names / snp2_gene_names)
        #    vs cluster loci col (gene-name loci on cluster-file side)
        elif loci_candidates & lookup['loci']:
            hit = 1
            method = 'loci_loci'

        # 3. loci → feature
        #    our gene-name tokens vs cluster feature col
        #    (cluster file lists a gene symbol rather than an rsID as the
        #     feature identifier — e.g. some loci-based cluster maps)
        elif loci_candidates & lookup['features']:
            hit = 1
            method = 'loci_feature'

        # 4. feature → loci
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

def main(pheno_data, filter_strategy='tiered',
         raw_feature_scores=None,
         use_set='holdout',
         balance_strategy='unbalanced',
         n_iterations=1000,
         alpha=0.05,
         min_effect_size=0.1,
         majority_threshold=0.6,
         delta_threshold=0.1,
         use_all=False):

    scores_path = os.path.join(pheno_data, 'scores')
    output_path = os.path.join(scores_path, 'functionalEnrichment')
    os.makedirs(output_path, exist_ok=True)

    # -------------------------------------------------------------------------
    # Load all important features from the filtered tiered CSV
    # (HighCases + LowControls + ExclusiveHighCases)
    # -------------------------------------------------------------------------
    print("\n--- Step 1: Loading filtered important features ---")
    filtered_df = load_genomic_features(scores_path, filter_strategy)
    if filtered_df.empty:
        print("[ERROR] No filtered features found. Ensure "
              "importantFeaturesAcrossCohortsAndTrainingData.Filtered."
              f"{filter_strategy}.csv exists.")
        return

    # -------------------------------------------------------------------------
    # Load feature annotation (for SNP expansion)
    # -------------------------------------------------------------------------
    print("\n--- Step 2: Loading feature annotation ---")
    annotation_df = load_feature_annotation(scores_path)

    # -------------------------------------------------------------------------
    # Clinical feature analysis (optional — only when raw_feature_scores given)
    # -------------------------------------------------------------------------
    clinical_df = pd.DataFrame()
    if raw_feature_scores:
        print(f"\n--- Step 2a: Clinical feature cross-cohort analysis ---")
        print(f"  Strategy: {balance_strategy}  |  alpha={alpha}  "
              f"|  min_effect={min_effect_size}")

        # Load PRS bins to define exclusive groups
        if use_set == 'holdout':
            prs_file = os.path.join(scores_path, 'combinedPRSGroups.holdout.filtered.csv')
            if not os.path.exists(prs_file):
                prs_file = os.path.join(scores_path, 'combinedPRSGroups.filtered.csv')
        else:
            prs_file = os.path.join(scores_path, 'combinedPRSGroups.filtered.csv')
        combined_prs = pd.read_csv(prs_file)
        bin_cols = [c for c in combined_prs.columns
                    if c.startswith('bin_') and 'combined' not in c]
        if not use_all:
            bin_cols = [c for c in bin_cols if 'all' not in c]
        combined_prs = combined_prs[['IID', 'PHENOTYPE'] + bin_cols].set_index('IID')
        cohort_groups = identify_exclusive_groups(combined_prs, bin_cols)

        # Load and normalise raw clinical features
        raw_df = pd.read_csv(raw_feature_scores)
        if 'Participant ID' in raw_df.columns:
            raw_df.rename(columns={'Participant ID': 'IID'}, inplace=True)
        raw_df.set_index('IID', inplace=True)
        cohort_names = [c.replace('bin_', '') for c in bin_cols]
        exclude_meta = {'PHENOTYPE', 'SEX', 'age'} | set(bin_cols) | set(cohort_names)
        feature_cols_raw = [c for c in raw_df.columns if c not in exclude_meta]

        canon_map, _ = normalise_series_names(feature_cols_raw)
        canonical_seen = {}
        rename_dict = {}
        for orig in feature_cols_raw:
            can = canon_map[orig]
            if can not in canonical_seen:
                canonical_seen[can] = orig
                rename_dict[orig] = can
            else:
                print(f"  [NORM] Dropping duplicate canonical '{can}': "
                      f"keeping '{canonical_seen[can]}', dropping '{orig}'")
        raw_df = raw_df[[c for c in raw_df.columns
                         if c in exclude_meta or c in rename_dict]]
        raw_df.rename(columns=rename_dict, inplace=True)
        feature_cols = [c for c in raw_df.columns if c not in exclude_meta]
        print(f"  Clinical columns: {len(feature_cols_raw)} raw "
              f"→ {len(feature_cols)} after normalisation")

        # Run MWU for each cohort
        cohort_results = {}
        for cohort, groups in cohort_groups.items():
            cases    = [i for i in groups['cases']    if i in raw_df.index]
            controls = [i for i in groups['controls'] if i in raw_df.index]
            print(f"  {cohort}: {len(cases)} cases | {len(controls)} controls")
            if len(cases) < 5 or len(controls) < 5:
                print(f"    [SKIP] too few samples")
                cohort_results[cohort] = pd.DataFrame()
                continue
            if balance_strategy == 'balanced':
                res = find_clinical_features_balanced(
                    raw_df, cases, controls, feature_cols,
                    alpha=alpha, min_effect_size=min_effect_size,
                    n_iterations=n_iterations,
                )
            else:
                res = find_clinical_features_unbalanced(
                    raw_df, cases, controls, feature_cols,
                    alpha=alpha, min_effect_size=min_effect_size,
                )
            cohort_results[cohort] = res

        # Apply cross-cohort specificity filter
        clinical_df = filter_clinical_features_cross_cohort(
            cohort_results,
            majority_threshold=majority_threshold,
            delta_threshold=delta_threshold,
        )
    else:
        print("\n  [skip] --raw_feature_scores not provided; "
              "clinical analysis omitted.")

    # -------------------------------------------------------------------------
    # Build per-cohort feature rows from the filtered CSV
    # -------------------------------------------------------------------------
    print("\n--- Step 3: Building per-cohort feature list ---")
    all_feature_rows = []

    for cohort in sorted(c for c in filtered_df['cohort'].unique() if pd.notna(c)):
        cohort_feats = filtered_df[filtered_df['cohort'] == cohort].copy()
        cohort_feats['raw_feature'] = cohort_feats['feature']
        cohort_feats['p_value']      = np.nan
        cohort_feats['p_value_fdr']  = np.nan
        cohort_feats['effect_size_r'] = np.nan
        n_genomic  = (cohort_feats['feature_source'] == 'genomic_shap').sum()
        n_clinical = (cohort_feats['feature_source'] == 'clinical').sum()
        print(f"  {cohort}: {n_genomic} genomic, {n_clinical} clinical features "
              f"(from filtered CSV)")
        all_feature_rows.append(cohort_feats)

    # Add cross-cohort clinical features (deduplicated against filtered CSV clinical rows)
    if not clinical_df.empty:
        # Avoid double-counting clinical features already present in filtered CSV
        existing_clinical = set(
            zip(filtered_df.loc[filtered_df['feature_source'] == 'clinical', 'feature'],
                filtered_df.loc[filtered_df['feature_source'] == 'clinical', 'cohort'])
        )
        novel_clinical = clinical_df[
            ~clinical_df.apply(
                lambda r: (r['feature'], r['cohort']) in existing_clinical, axis=1
            )
        ].copy()
        print(f"\n  Clinical cross-cohort: {len(clinical_df)} features "
              f"→ {len(novel_clinical)} not already in filtered CSV")
        if not novel_clinical.empty:
            all_feature_rows.append(novel_clinical)

    if not all_feature_rows:
        print("\n[ERROR] No important features found. Check inputs and thresholds.")
        return

    # -------------------------------------------------------------------------
    # Combine, deduplicate, and expand to one row per SNP
    # -------------------------------------------------------------------------
    keep_cols = ['raw_feature', 'feature', 'cohort', 'feature_source',
                 'p_value', 'p_value_fdr', 'effect_size_r',
                 'odds_ratio', 'specificity_tier',
                 'training_data_used', 'n_exclusive_cohorts', 'exclusive_cohorts']

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

    # ── Issue 2: re-infer feature_source after expansion ─────────────────────
    # A raw_feature like "Total bilirubin,rs459193,rs75981964" is classified as
    # 'clinical' by _infer_feature_source (because not all comma-parts are rsIDs).
    # After expansion each part becomes its own row: the rsID rows must be
    # reclassified as 'genomic_shap' now that their resolved feature is a bare rsID.
    expanded_df['feature_source'] = expanded_df['feature'].apply(_infer_feature_source)

    # Deduplicate expanded (same snp can appear in snp1 of one row and snp2 of another)
    expanded_df.drop_duplicates(
        subset=['feature', 'cohort', 'feature_source'], inplace=True
    )
    print(f"Features after SNP expansion: {len(expanded_df)}")

    # Gene-name coverage diagnostic — loci matching only works when gene names are present
    genomic_rows = expanded_df[expanded_df['feature_source'] == 'genomic_shap']
    has_genes = genomic_rows['snp1_gene_names'].notna() | genomic_rows['snp2_gene_names'].notna()
    print(f"  Gene-name coverage: {has_genes.sum()}/{len(genomic_rows)} genomic rows "
          f"have ≥1 gene name from feature_scores_annotated.csv")

    # -------------------------------------------------------------------------
    # Load functional cluster files
    # -------------------------------------------------------------------------
    print("\n--- Step 4: Loading functional cluster files ---")
    # NOTE: Mahajan map is applied in Step 4b (after subtype_path is set).
    subtype_path = os.path.join(os.path.dirname(pheno_data), 'subtypeResearch')

    # Populate clinical alias map BEFORE any normalise_name calls on cluster data
    global _CLINICAL_FEATURE_MAP
    _CLINICAL_FEATURE_MAP = load_clinical_feature_map(subtype_path)

    # Load Mahajan SNP → clinical-feature bridge map
    mahajan_map = load_mahajan_map(subtype_path)

    # Load ALL CSVs in subtypeResearch/ as potential cluster maps, excluding
    # the two files that are loaded separately for other purposes.
    # A pattern like *functional_cluster*.csv misses files named differently,
    # e.g. kim_function_cluster_*.csv (no 'al') and multerer_fc_*.csv ('fc'
    # abbreviation), so we enumerate all CSVs and explicitly exclude the known
    # non-cluster files instead.
    _EXCLUDE_FROM_CLUSTER = {
        'clinical_feature_map.csv',        # loaded by load_clinical_feature_map()
        'mahajan_clinical_loci_map.csv',   # loaded by load_mahajan_map()
        'multerer_fc_clinical_trait_map.csv',  # old results, not used
        'multerer_fc_loci_map.csv',            # old results, not used
    }
    cluster_files = sorted(
        f for f in glob.glob(os.path.join(subtype_path, '*.csv'))
        if os.path.basename(f) not in _EXCLUDE_FROM_CLUSTER
    )

    if not cluster_files:
        print(f"  [WARN] No cluster files found in {subtype_path}")

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
    # Step 4b: Annotate genomic features with Mahajan clinical bridge
    # -------------------------------------------------------------------------
    print("\n--- Step 4b: Applying Mahajan SNP→clinical-feature linkage ---")
    expanded_df = apply_mahajan_clinical_features(expanded_df, mahajan_map)
    n_mahajan = (expanded_df['mahajan_clinical_feature'] != '').sum()
    print(f"  Genomic rows with a Mahajan clinical link: {n_mahajan}/{len(expanded_df)}")

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
            mahajan_clinical_feature=row.get('mahajan_clinical_feature', ''),
        )
        # Count how many functional cluster columns have a positive hit (1)
        n_matched = sum(cluster_hits.values())

        output_rows.append({
            'raw_feature':               row.get('raw_feature', feat_val),
            'feature':                   feat_val,
            'cohort':                    row['cohort'],
            'training_data_used':        row.get('training_data_used', np.nan),
            'feature_source':            row.get('feature_source'),
            'p_value':                   row.get('p_value'),
            'p_value_fdr':               row.get('p_value_fdr'),
            'effect_size_r':             row.get('effect_size_r'),
            'odds_ratio':                row.get('odds_ratio'),
            'specificity_tier':          row.get('specificity_tier'),
            'n_exclusive_cohorts':       row.get('n_exclusive_cohorts', np.nan),
            'exclusive_cohorts':         row.get('exclusive_cohorts', np.nan),
            'snp1_gene_names':           row.get('snp1_gene_names'),
            'snp2_gene_names':           row.get('snp2_gene_names'),
            'mahajan_clinical_feature':  row.get('mahajan_clinical_feature', ''),
            'mahajan_source':            row.get('mahajan_source', ''),
            'match_method':              match_method,
            'n_clusters_matched':        n_matched,
            **cluster_hits,              # binary 0/1 per functional cluster column
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

    # -------------------------------------------------------------------------
    # HLA autoimmune cluster: features whose name starts with HLA- are assigned
    # directly to the hla_autoimmune cluster (→ SAID in cluster_subtype_map.csv).
    # This bypasses the SNP/locus-matching pipeline because HLA allele calls
    # (e.g. HLA-drb1*401) are not present in the standard functional cluster CSVs.
    # -------------------------------------------------------------------------
    HLA_COL = 'hla_cluster__HLA Region'
    hla_mask = output_df['feature'].str.match(r'^HLA[-_*]', case=False, na=False)
    output_df[HLA_COL] = hla_mask.astype(int)
    if hla_mask.any():
        print(f"  HLA autoimmune cluster ({HLA_COL}): {hla_mask.sum()} features assigned")
        # Update match_method for HLA rows that had no other match
        if 'match_method' in output_df.columns:
            output_df.loc[
                hla_mask & (output_df['match_method'] == 'none'),
                'match_method'
            ] = 'hla_direct'
    all_cluster_col_keys = all_cluster_col_keys + [HLA_COL]

    # Recompute n_clusters_matched after enforcing binary ints
    if all_cluster_col_keys:
        output_df['n_clusters_matched'] = output_df[all_cluster_col_keys].sum(axis=1)

    # Canonical column order: metadata | match tracking | cluster binary columns
    meta_cols = [
        'raw_feature', 'feature', 'cohort', 'training_data_used', 'feature_source',
        'p_value', 'p_value_fdr', 'effect_size_r', 'odds_ratio',
        'specificity_tier',
        # ExclusiveHighCases cross-cohort prevalence annotations
        'n_exclusive_cohorts', 'exclusive_cohorts',
        'snp1_gene_names', 'snp2_gene_names',
        # Mahajan clinical linkage (populated for genomic features only)
        'mahajan_clinical_feature', 'mahajan_source',
        'match_method', 'n_clusters_matched',
    ]
    # Only keep meta_cols that were actually populated (some are ExclusiveHighCases-only)
    present_meta = [c for c in meta_cols if c in output_df.columns]
    cluster_cols_out = [c for c in output_df.columns if c not in meta_cols]
    output_df = output_df[present_meta + cluster_cols_out]

    # -------------------------------------------------------------------------
    # Save and summarise
    # -------------------------------------------------------------------------
    out_file = os.path.join(
        output_path,
        f'functional_cluster_enrichment_{filter_strategy}.csv'
    )
    output_df.to_csv(out_file, index=False)

    n_any_cluster = int((output_df['n_clusters_matched'] > 0).sum()) \
        if 'n_clusters_matched' in output_df.columns else 0
    print(f"\nSaved → {out_file}")
    print(f"Shape : {output_df.shape}")
    print(f"Features mapped to ≥1 cluster: {n_any_cluster}/{len(output_df)}")

    print("\nPer-cohort feature counts (training_data_used × cohort × source):")
    grp_cols = [c for c in ['training_data_used', 'cohort', 'feature_source']
                if c in output_df.columns]
    summary = (
        output_df.groupby(grp_cols)
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
        '--filter_strategy', default='tiered',
        choices=['strict', 'majority', 'differential', 'tiered'],
        help='Filtering strategy used in importantFeaturesAcrossCohortsAndTrainingData.Filtered.*.csv (default: tiered)',
    )
    parser.add_argument(
        '--raw_feature_scores',
        help='Path to raw (untransformed) clinical features CSV. '
             'When provided, clinical features are tested with Mann-Whitney U '
             'across all cohorts and filtered for cross-cohort specificity.',
    )
    parser.add_argument(
        '--use_set', default='holdout',
        choices=['holdout', 'validation', 'all'],
        help='PRS bin file to use for defining exclusive case/control groups '
             'in clinical analysis (default: holdout)',
    )
    parser.add_argument(
        '--balance_strategy', default='unbalanced',
        choices=['unbalanced', 'balanced'],
        help='unbalanced (default): all controls vs exclusive cases. '
             'balanced: bootstrap n_cases controls over n_iterations.',
    )
    parser.add_argument(
        '--n_iterations', type=int, default=1000,
        help='Bootstrap iterations for balanced clinical strategy (default: 1000)',
    )
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='FDR significance threshold for clinical features (default: 0.05)',
    )
    parser.add_argument(
        '--min_effect_size', type=float, default=0.1,
        help='Minimum |rank-biserial r| for clinical features (default: 0.1)',
    )
    parser.add_argument(
        '--majority_threshold', type=float, default=0.6,
        help='Fraction of other cohorts that must be non-significant for '
             'tier-2 clinical specificity (default: 0.6)',
    )
    parser.add_argument(
        '--delta_threshold', type=float, default=0.1,
        help='Minimum |r| advantage over mean of other cohorts for '
             'tier-3 clinical specificity (default: 0.1)',
    )
    parser.add_argument(
        '--use_all', action='store_true',
        help='Include "all" cohort bins in clinical analysis (excluded by default)',
    )

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get('PHENO_DATA')
    raw_feature_scores = args.raw_feature_scores or os.environ.get('RAW_FEATURE_SCORES')

    if not pheno_data:
        parser.error('Provide --pheno_data or set PHENO_DATA env var')

    print(f"pheno_data        : {pheno_data}")
    print(f"filter_strategy   : {args.filter_strategy}")
    print(f"raw_feature_scores : {raw_feature_scores or '(not provided — clinical analysis skipped)'}")
    if raw_feature_scores:
        print(f"use_set           : {args.use_set}")
        print(f"balance_strategy  : {args.balance_strategy}")
        print(f"n_iterations      : {args.n_iterations}")
        print(f"alpha             : {args.alpha}")
        print(f"min_effect_size   : {args.min_effect_size}")
        print(f"majority_threshold: {args.majority_threshold}")
        print(f"delta_threshold   : {args.delta_threshold}")
        print(f"use_all           : {args.use_all}")

    main(
        pheno_data=pheno_data,
        filter_strategy=args.filter_strategy,
        raw_feature_scores=raw_feature_scores,
        use_set=args.use_set,
        balance_strategy=args.balance_strategy,
        n_iterations=args.n_iterations,
        alpha=args.alpha,
        min_effect_size=args.min_effect_size,
        majority_threshold=args.majority_threshold,
        delta_threshold=args.delta_threshold,
        use_all=args.use_all,
    )
