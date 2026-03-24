#!/usr/bin/env python3
"""
prs_bnmf_subtype_analysis.py

Performs Bayesian-inspired Non-negative Matrix Factorization (bNMF) on ePRS
model scores combined with clinical measures to identify phenotypic subtypes
within disease risk cohorts.

Methodology inspired by Udler et al. (2018) PLOS Medicine which identified
5 T2D subtypes using soft clustering of genetic risk scores and clinical
biomarkers. DOI: https://doi.org/10.1371/journal.pmed.1002654

Two input modes
---------------
  prs_scores        Uses aggregate scaled_prs value from each *.mixed.prs.csv
                    model combined with clinical measures (analogous to Udler's
                    use of pathway-level genetic risk scores).

  weighted_features Uses the per-individual weighted SNP/interaction
                    contributions from *.mixed.prs.csv files (the raw columns
                    before summing to prs). Only the important features listed
                    in --important_features_file are retained to keep the
                    matrix tractable.

Outputs
-------
  {pheno_data}/scores/bnmf/
      cluster_assignments.csv       soft + hard cluster membership per IID
      feature_loadings.csv          H matrix (feature weights per cluster)
      basis_matrix.csv              W matrix (individual cluster memberships)
      prs_clinical_associations.csv Spearman correlations, FDR-corrected
      cophenetic_curve.csv          cophenetic correlation per k tested

  {pheno_data}/figures/bnmf/
      cophenetic_curve.png
      consensus_map_k{k}.png        (for each k tested)
      cluster_profiles.png
      prs_clinical_heatmap.png
      cluster_udler_comparison.png  (T2D validation only)
"""

import os
import sys
import gc
import glob
import re
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, cophenet, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler, QuantileTransformer
from sklearn.impute import SimpleImputer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns

warnings.filterwarnings("ignore")


# ============================================================================
# CONSTANTS
# ============================================================================

# Udler et al. 2018 T2D subtype reference profiles.
# Each entry maps standardized clinical feature keywords to expected direction
# relative to population mean (+1 higher, -1 lower, 0 neutral).
# Used only for T2D phenotype validation.
UDLER_PROFILES = {
    'SIDD': {   # Severe Insulin-Deficient Diabetes
        'hba1c': +1, 'bmi': -1, 'c_peptide': -1, 'hdl': 0,
        'triglyceride': 0, 'age': -1,
        'description': 'Low insulin secretion, low BMI, high HbA1c'
    },
    'SIRD': {   # Severe Insulin-Resistant Diabetes
        'hba1c': +1, 'bmi': +1, 'c_peptide': +1, 'hdl': -1,
        'triglyceride': +1, 'age': 0,
        'description': 'High insulin resistance, high BMI, dyslipidemia'
    },
    'MOD': {    # Mild Obesity-related Diabetes
        'hba1c': 0, 'bmi': +1, 'c_peptide': 0, 'hdl': 0,
        'triglyceride': 0, 'age': 0,
        'description': 'Mild hyperglycemia, obesity-driven'
    },
    'MARD': {   # Mild Age-Related Diabetes
        'hba1c': 0, 'bmi': 0, 'c_peptide': 0, 'hdl': 0,
        'triglyceride': 0, 'age': +1,
        'description': 'Mild diabetes, older age at onset'
    },
    'SAID': {   # Severe Autoimmune Diabetes
        'hba1c': +1, 'bmi': -1, 'c_peptide': -1, 'hdl': 0,
        'triglyceride': 0, 'age': -1,
        'description': 'Autoimmune-mediated, low BMI, high HbA1c'
    },
}

# Keyword mapping from raw column names to standardized clinical concepts
# (lower-case substring matching)
CLINICAL_KEYWORDS = {
    'bmi': 'bmi',
    'body mass': 'bmi',
    'hba1c': 'hba1c',
    'haemoglobin a1c': 'hba1c',
    'glycated haemoglobin': 'hba1c',
    'glycohaemoglobin': 'hba1c',
    'c-peptide': 'c_peptide',
    'cpeptide': 'c_peptide',
    'hdl': 'hdl',
    'triglyceride': 'triglyceride',
    'age at recruit': 'age',
    'age at diagnosis': 'age',
    'systolic': 'systolic_bp',
    'diastolic': 'diastolic_bp',
    'glucose': 'glucose',
    'insulin': 'insulin',
    'cholesterol': 'cholesterol',
    'creatinine': 'creatinine',
    'egfr': 'egfr',
    'waist': 'waist',
    'alt': 'alt',
    'ast': 'ast',
    'ldl': 'ldl',
}

# Maximum number of individuals used to build the consensus matrix for
# visualisation only.  The full n×n matrix is O(n²) memory; at n=67k that is
# ~18 GB per k even in float32, which is unplottable anyway.  A random
# subsample of this size gives a clear, representative heatmap.
MAX_CONSENSUS_PLOT_N = 5000

# Okabe-Ito colorblind-friendly palette
CLUSTER_PALETTE = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7', '#999999',
    '#44AA99', '#332288',
]

UDLER_COLORS = {
    'SIDD': '#66C2A5',
    'SIRD': '#FC8D62',
    'MOD': '#8DA0CB',
    'MARD': '#E78AC3',
    'SAID': '#A6D854',
    'unassigned': '#CCCCCC',
}


# ============================================================================
# DATA LOADING
# ============================================================================

def _read_prs_file(filepath):
    """Read a single *.mixed.prs.csv and return a DataFrame indexed by IID."""
    df = pd.read_csv(filepath)
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'IID'}, inplace=True)
    if 'IID' not in df.columns:
        raise ValueError(f"No IID column found in {filepath}")
    df = df.set_index('IID')
    return df


def _model_name_from_path(filepath, models_to_use):
    """Extract model name from a PRS file path using longest-match regex."""
    filename = os.path.basename(filepath)
    models_sorted = sorted(models_to_use, key=len, reverse=True)
    for model in models_sorted:
      pattern = rf'^{re.escape(model)}(\..+)?\.mixed\.prs\.csv$'
      if re.search(pattern, filename):
            return model
    return None


def load_prs_files(scores_path, models_to_use=None):
    """
    Discover and load *.mixed.prs.csv files (excluding holdout).

    Returns
    -------
    dict  { model_name: pd.DataFrame }
        DataFrames indexed by IID. Columns are per-feature weights plus
        prs, scaled_prs, high_risk, PHENOTYPE.
    """
    pattern = os.path.join(scores_path, 'prsScores', '*mixed*.csv')
    all_files = glob.glob(pattern)
    prs_files = [f for f in all_files if 'holdout' not in os.path.basename(f).lower()]

    if models_to_use:
        prs_files = [f for f in prs_files
                     if _model_name_from_path(f, models_to_use) is not None]

    prs_dict = {}
    for filepath in prs_files:
        if models_to_use:
            name = _model_name_from_path(filepath, models_to_use)
        else:
            # Derive from filename stem before first dot
            name = os.path.basename(filepath).split('.')[0]
        if name is None:
            continue
        try:
            prs_dict[name] = _read_prs_file(filepath)
            print(f"  Loaded PRS file: {os.path.basename(filepath)} → model '{name}'")
        except Exception as exc:
            print(f"  WARNING: could not load {filepath}: {exc}")

    return prs_dict


def build_prs_score_matrix(prs_dict):
    """
    Mode: prs_scores
    Extract only the scaled_prs column per model and merge into a single
    wide matrix (analogous to Udler's pathway-level genetic risk scores).

    Returns
    -------
    pd.DataFrame  columns = {model}_prs, index = IID
    pd.Series     PHENOTYPE per IID
    """
    frames = []
    phenotype_series = None

    for model_name, df in prs_dict.items():
        col_name = f'{model_name}_prs'
        if 'scaled_prs' in df.columns:
            frames.append(df[['scaled_prs']].rename(columns={'scaled_prs': col_name}))
        elif 'prs' in df.columns:
            frames.append(df[['prs']].rename(columns={'prs': col_name}))
        else:
            print(f"  WARNING: no prs column for model '{model_name}', skipping.")

        if phenotype_series is None and 'PHENOTYPE' in df.columns:
            phenotype_series = df['PHENOTYPE']

    if not frames:
        raise ValueError("No PRS score columns could be extracted.")

    merged = frames[0]
    for f in frames[1:]:
        merged = merged.join(f, how='outer')

    return merged, phenotype_series


def build_weighted_feature_matrix(prs_dict, important_features=None):
    """
    Mode: weighted_features
    Concatenate per-individual weighted feature columns across all PRS models.
    Feature columns are the raw SNP/interaction contribution values
    (individual dosage × coefficient) before summation to prs.

    If important_features (pd.DataFrame with 'feature' and 'model' columns)
    is provided, restrict to those features only.

    Returns
    -------
    pd.DataFrame  columns = {model}__{feature}, index = IID
    pd.Series     PHENOTYPE per IID
    """
    frames = []
    phenotype_series = None
    prs_cols = {'prs', 'scaled_prs', 'high_risk', 'PHENOTYPE', 'centile_bin'}

    for model_name, df in prs_dict.items():
        feature_cols = [c for c in df.columns
                        if c not in prs_cols and not c.startswith('PC')]

        if important_features is not None and 'feature' in important_features.columns:
            model_feats = important_features[
                important_features['model'] == model_name
            ]['feature'].tolist()
            feature_cols = [c for c in feature_cols if c in model_feats]

        if not feature_cols:
            print(f"  WARNING: no feature columns for model '{model_name}' after filtering.")
            continue

        sub = df[feature_cols].copy()
        sub.columns = [f'{model_name}__{c}' for c in feature_cols]
        frames.append(sub)

        if phenotype_series is None and 'PHENOTYPE' in df.columns:
            phenotype_series = df['PHENOTYPE']

    if not frames:
        raise ValueError("No weighted feature columns could be extracted.")

    merged = frames[0]
    for f in frames[1:]:
        merged = merged.join(f, how='outer')

    return merged, phenotype_series


def load_clinical_measures(clinical_file, iid_col='eid'):
  """
  Load participant_environment.csv (or any clinical measures CSV).

  Returns
  -------
  pd.DataFrame  indexed by IID
  """
  df = pd.read_csv(clinical_file, low_memory=False)
  
  # Drop duplicate columns (keep first occurrence)
  n_dupes = df.columns.duplicated().sum()
  if n_dupes:
    print(f"  Warning: dropping {n_dupes} duplicate column(s)")
    df = df.loc[:, ~df.columns.duplicated()]
    
  id_candidates = [iid_col, 'IID', 'Participant ID', 'participant_id', 'f.eid']
  id_col_found = None
  for c in id_candidates:
    if c in df.columns:
      id_col_found = c
      break
  if id_col_found is None:
    raise ValueError(
      f"Cannot identify ID column in {clinical_file}. "
      f"Expected one of: {id_candidates}"
    )
  df = df.set_index(id_col_found)
  df.index.name = 'IID'
  
  # Keep only numeric columns
  numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
  print(f"  Loaded clinical measures: {df.shape[0]} individuals, "
      f"{len(numeric_cols)} numeric features")
  return df[numeric_cols]


def map_feature_to_concept(col_name):
    """Map a raw column name to a standardized clinical concept string."""
    lower = col_name.lower()
    for keyword, concept in CLINICAL_KEYWORDS.items():
        if keyword in lower:
            return concept
    return None


# ============================================================================
# PREPROCESSING
# ============================================================================

def impute_and_scale(df, imputer_strategy='median', n_neighbors=5,
                     scale_method='quantile'):
    """
    Impute missing values (median) then scale to [0, 1] for NMF.
    Columns with >80% missing are dropped before imputation.

    Column-type-aware scaling (detected from column-name prefix):

      gen_* (genomic features — dosage × beta coefficient)
          Per-column QuantileTransformer (uniform [0,1] by rank).
          Each SNP contribution is ranked within the subgroup so the most
          risk-increasing genotype maps to 1.0 and the most risk-reducing
          (including negative/protective contributions) maps to 0.0.
          This eliminates the sparsity collapse caused by clipping negatives
          to 0: with p99-clipping, protective SNPs (negative weight) have
          zero variance after clipping and are filtered out; positive-weight
          SNPs with low MAF have 80%+ zeros and also collapse NMF.

      clin_* and all other columns (clinical measures, PRS anchors)
          scale_method='quantile' (default): QuantileTransformer maps each
          feature to uniform [0,1] by rank, removing the shared "average
          individual" direction that dominates MinMaxScaler in homogeneous
          subgroups and causes NMF cluster collapse.
          scale_method='minmax': MinMaxScaler; preserves raw distributions
          but amplifies correlated features.

    Returns
    -------
    pd.DataFrame  scaled, no NaNs, same index as input
    list          columns retained after missingness filter
    """
    miss_frac = df.isnull().mean()
    keep_cols = miss_frac[miss_frac <= 0.8].index.tolist()
    dropped = miss_frac[miss_frac > 0.8].index.tolist()
    if dropped:
        print(f"  Dropped {len(dropped)} columns with >80% missing: "
              f"{dropped[:5]}{'...' if len(dropped) > 5 else ''}")
    df = df[keep_cols]

    imputer = SimpleImputer(strategy=imputer_strategy)
    arr_imp = imputer.fit_transform(df.values)
    df_imp  = pd.DataFrame(arr_imp, index=df.index, columns=keep_cols)

    # ── Split columns by type ───────────────────────────────────────────────
    gen_cols   = [c for c in keep_cols if c.startswith('gen_')]
    other_cols = [c for c in keep_cols if not c.startswith('gen_')]

    result = pd.DataFrame(index=df.index, columns=keep_cols, dtype=float)

    # ── Genomic: per-column quantile normalisation ─────────────────────────
    # Rank each SNP contribution within the subgroup → uniform [0, 1].
    # Negatives (protective variants) are included in the ranking rather than
    # clipped to 0; this removes the sparsity that causes NMF to collapse.
    if gen_cols:
        n_q = min(1000, df_imp.shape[0])
        qt_gen = QuantileTransformer(
            output_distribution='uniform',
            n_quantiles=n_q,
            random_state=42,
        )
        result[gen_cols] = qt_gen.fit_transform(df_imp[gen_cols].values)
        neg_frac = (df_imp[gen_cols].values < 0).mean()
        print(f"    Genomic scaling: per-column quantile normalisation "
              f"({len(gen_cols)} gen_ columns; "
              f"{neg_frac*100:.1f}% of raw values were negative — "
              f"now included in ranking)")

    # ── Clinical / PRS anchors: per-column quantile or minmax ──────────────
    if other_cols:
        other_arr = df_imp[other_cols].values
        if scale_method == 'quantile':
            n_q = min(1000, other_arr.shape[0])
            scaler = QuantileTransformer(
                output_distribution='uniform',
                n_quantiles=n_q,
                random_state=42,
            )
        else:
            scaler = MinMaxScaler()
        result[other_cols] = scaler.fit_transform(other_arr)

    return result, keep_cols


# ============================================================================
# ASSOCIATION ANALYSIS (PRS models ↔ clinical measures)
# ============================================================================

def compute_prs_clinical_associations(prs_matrix, clinical_matrix):
    """
    Compute Spearman correlation between each PRS model column and each
    clinical measure. FDR correction (Benjamini-Hochberg) across all tests.

    Returns
    -------
    pd.DataFrame  rows = clinical features, cols = PRS models,
                  values = rho (Spearman r), with p-value columns
    pd.DataFrame  full association table with FDR-corrected q-values
    """
    shared_idx = prs_matrix.index.intersection(clinical_matrix.index)
    print(f"  Association analysis: {len(shared_idx)} shared individuals")

    prs_sub = prs_matrix.loc[shared_idx]
    clin_sub = clinical_matrix.loc[shared_idx]

    records = []
    for prs_col in prs_sub.columns:
        for clin_col in clin_sub.columns:
            x = prs_sub[prs_col]
            y = clin_sub[clin_col]
            mask = x.notna() & y.notna()
            if mask.sum() < 30:
                continue
            rho, pval = spearmanr(x[mask], y[mask])
            records.append({
                'prs_model': prs_col,
                'clinical_feature': clin_col,
                'rho': rho,
                'p_value': pval,
                'n': int(mask.sum()),
            })

    assoc_df = pd.DataFrame(records)
    if assoc_df.empty:
        print("  WARNING: no associations computed (too few shared individuals?).")
        return pd.DataFrame(), pd.DataFrame()

    _, q_values, _, _ = multipletests(assoc_df['p_value'], method='fdr_bh')
    assoc_df['q_value'] = q_values
    assoc_df['significant'] = assoc_df['q_value'] < 0.05

    # Pivot to rho matrix for heatmap
    rho_pivot = assoc_df.pivot(
        index='clinical_feature', columns='prs_model', values='rho'
    )
    return rho_pivot, assoc_df


def plot_association_heatmap(rho_pivot, output_path, title='PRS–Clinical Spearman ρ'):
    """Save a clustered heatmap of PRS × clinical Spearman correlations."""
    if rho_pivot.empty:
        return
    fig_h = max(6, len(rho_pivot) * 0.3)
    fig_w = max(6, len(rho_pivot.columns) * 0.8)
    plt.figure(figsize=(fig_w, fig_h))
    sns.clustermap(
        rho_pivot.fillna(0),
        cmap='RdBu_r', center=0, vmin=-0.5, vmax=0.5,
        linewidths=0.3, linecolor='#cccccc',
        annot=len(rho_pivot) <= 30,
        fmt='.2f' if len(rho_pivot) <= 30 else '',
        figsize=(fig_w, fig_h),
    )
    plt.suptitle(title, y=1.01, fontsize=11)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved association heatmap → {output_path}")


# ============================================================================
# bNMF CONSENSUS CLUSTERING
# ============================================================================

def _run_single_nmf(X, k, random_state, l1_ratio=0.1, alpha_W=0.1, alpha_H=0.1,
                    max_iter=500):
    """
    Single NMF run with KL divergence (Kullback-Leibler) + sparsity
    regularization on both W and H.

    alpha_H must be > 0 when features include sparse/binary-like columns
    (e.g., HLA allele dosage contributions).  Without regularisation on H,
    the MU solver can compensate for L1-shrunk W by growing H arbitrarily
    large, producing degenerate components where almost all individuals are
    assigned to cluster_1 (the catch-all component).

    Returns W (n_samples × k) and H (k × n_features).
    """
    model = NMF(
        n_components=k,
        init='random',          # random init allows 30 runs to explore different
                                # solutions; nndsvda is deterministic-near-SVD and
                                # locks all runs near the rank-1 fixed point for
                                # homogeneous subgroups (root cause of collapse)
        beta_loss='kullback-leibler',
        solver='mu',
        max_iter=max_iter,
        random_state=random_state,
        l1_ratio=l1_ratio,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
        tol=1e-4,
    )
    W = model.fit_transform(X).astype(np.float32)
    H = model.components_.astype(np.float32)
    return W, H, model.reconstruction_err_


def _soft_consensus_matrix(W_list):
    """
    Build consensus matrix from a list of W matrices (soft memberships).
    Entry C[i,j] = mean cosine similarity of individual membership vectors.
    Averages over all NMF runs to give a stable pairwise co-clustering measure.
    """
    n = W_list[0].shape[0]
    C = np.zeros((n, n), dtype=np.float32)
    for W in W_list:
        # Normalize rows to unit L1 (soft membership fractions)
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        W_norm = W / row_sums
        # Outer product: C[i,j] = dot(w_i, w_j) — measures co-cluster affinity
        C += W_norm @ W_norm.T
    C /= len(W_list)
    return C


def _cophenetic_correlation(C):
    """
    Compute cophenetic correlation coefficient of the consensus matrix C.
    Higher values indicate more stable clustering.

    Returns NaN when the consensus matrix is degenerate (near-constant), which
    happens when all NMF runs collapse to a single dominant component.  The
    caller should treat NaN as "unresolvable clustering for this k".
    """
    dist = 1 - C
    np.fill_diagonal(dist, 0)

    # Degenerate check: if all distances are effectively zero every individual
    # was always co-clustered → correlation is undefined.
    if dist.max() < 1e-9:
        return float('nan')

    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method='average')
    coph_dist = cophenet(Z, condensed)
    val = float(coph_dist[0])
    # scipy returns nan when condensed is constant (near-degenerate case)
    return val


def _pca_feature_selection(X_df, n_pcs=20, top_per_pc=5):
    """
    Select features proportional to each PC's explained variance ratio.

    Instead of a fixed top_per_pc allocation for every PC, the number of
    features drawn from PC_i is weighted by its share of the total explained
    variance across the n_pcs components:

        n_i = max(1, round(evr[i] / sum(evr) * total_budget))

    where total_budget = n_pcs * top_per_pc.

    Example with top_per_pc=5, n_pcs=10, total_budget=50:
        PC1 explains 20% of the n_pcs variance  →  ~10 features
        PC5 explains  8%                         →   ~4 features
        PC10 explains 2%                         →   ~1 feature

    High-variance PCs represent stronger biological signals and contribute
    more features; low-variance PCs contribute fewer but are still
    represented so rare axes of variation are not discarded entirely.

    Returns X_df restricted to the deduplicated union of selected features,
    preserving original column order.
    """
    from sklearn.decomposition import TruncatedSVD
    X = X_df.values.astype(np.float32)
    X_centered = X - X.mean(axis=0)
    n_comp = min(n_pcs, X.shape[1] - 1, X.shape[0] - 1)
    if n_comp < 1:
        return X_df

    svd = TruncatedSVD(n_components=n_comp, random_state=42)
    svd.fit(X_centered)

    evr          = svd.explained_variance_ratio_   # shape (n_comp,)
    total_evr    = evr.sum()
    total_budget = n_comp * top_per_pc

    selected  = set()
    alloc_log = []
    for i in range(n_comp):
        n_for_pc = max(1, round(evr[i] / total_evr * total_budget))
        top_idx  = np.argsort(np.abs(svd.components_[i]))[-n_for_pc:]
        selected.update(top_idx.tolist())
        alloc_log.append((i + 1, evr[i] * 100, n_for_pc))

    # Preserve original column order
    kept = [X_df.columns[i] for i in sorted(selected)]

    # Print per-PC allocation summary (top 5) so the user can audit
    summary = '  '.join(f'PC{pc}:{pct:.1f}%→{n}f'
                        for pc, pct, n in alloc_log[:5])
    if n_comp > 5:
        summary += f'  … (+{n_comp - 5} more PCs)'
    print(f"    PC feature selection: {len(kept)} features selected "
          f"(budget={total_budget}, proportional to variance)")
    print(f"      {summary}")
    return X_df[kept]


def run_consensus_bnmf(X, k_min=2, k_max=8, n_runs=30,
                       l1_ratio=0.1, alpha_W=0.1, alpha_H=0.1, max_iter=500,
                       verbose=True,
                       feature_select_method='variance',
                       n_pcs_for_features=20,
                       top_features_per_pc=5):
    """
    Consensus bNMF: for each k in [k_min, k_max], run NMF n_runs times,
    build the soft consensus matrix, and compute cophenetic correlation.

    Returns
    -------
    k_results : dict  {k: {'C': consensus_matrix, 'coph': float,
                            'W_list': [...], 'H_list': [...]}}
    coph_df   : pd.DataFrame  with columns k, cophenetic_correlation, mean_error
    """
    k_results = {}
    coph_records = []

    for k in range(k_min, k_max + 1):
        if verbose:
            print(f"  Running {n_runs} NMF restarts for k={k} ...", flush=True)
        W_list, H_list, errors = [], [], []

        for run in range(n_runs):
            try:
                W, H, err = _run_single_nmf(
                    X, k,
                    random_state=run * 17 + k,
                    l1_ratio=l1_ratio,
                    alpha_W=alpha_W,
                    alpha_H=alpha_H,
                    max_iter=max_iter,
                )
                W_list.append(W)
                H_list.append(H)
                errors.append(err)
            except Exception as exc:
                if verbose:
                    print(f"    run {run} failed: {exc}")

        if not W_list:
            print(f"  All runs failed for k={k}. Skipping.")
            continue

        # ── Collapse detection ──────────────────────────────────────────────
        # For each run, compute the fraction of individuals whose argmax W
        # column is the dominant (most-used) component.  Values near 1.0
        # mean every run put almost everyone in a single cluster.
        collapse_fracs = []
        for W in W_list:
            row_sums = W.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1
            W_norm = W / row_sums
            labels = np.argmax(W_norm, axis=1)
            dominant_frac = np.bincount(labels, minlength=k).max() / len(labels)
            collapse_fracs.append(dominant_frac)
        mean_collapse = float(np.mean(collapse_fracs))
        collapsed = mean_collapse > 0.90
        if collapsed and verbose:
            print(f"    *** k={k}: NMF COLLAPSE detected — {mean_collapse:.1%} of "
                  f"individuals assigned to dominant component on average across runs. "
                  f"Cophenetic correlation will be undefined (NaN). "
                  f"Consider reducing k_max, tightening feature filters, or "
                  f"increasing alpha_H.")

        # Build C only to compute cophenetic, then immediately free it.
        # C is O(n_samples²) and is the primary memory bottleneck.
        C = _soft_consensus_matrix(W_list)
        coph = _cophenetic_correlation(C)
        del C
        gc.collect()

        # Store W/H lists (O(n_samples × k) per run — much smaller than C).
        # C is NOT stored here; it will be recomputed one k at a time for
        # plotting after optimal k selection (see main).
        k_results[k] = {
            'coph': coph,
            'W_list': W_list,
            'H_list': H_list,
            'mean_collapse_frac': mean_collapse,
        }
        coph_records.append({
            'k': k,
            'cophenetic_correlation': coph,
            'mean_reconstruction_error': float(np.mean(errors)),
            'n_runs': len(W_list),
            'mean_collapse_frac': mean_collapse,
        })

        if verbose:
            coph_str = f"{coph:.4f}" if not np.isnan(coph) else "NaN (collapsed)"
            print(f"    k={k}: cophenetic={coph_str}  "
                  f"mean_error={np.mean(errors):.4f}  "
                  f"collapse={mean_collapse:.1%}")

    coph_df = pd.DataFrame(coph_records)
    return k_results, coph_df


def select_optimal_k(coph_df, method='elbow'):
    """
    Select optimal k from cophenetic correlation curve.

    method='max'   → k with highest cophenetic correlation
    method='elbow' → k at the "elbow" (largest second derivative drop)

    NaN cophenetic values (produced by degenerate / collapsed NMF runs) are
    excluded before selection.  If all k values are NaN the smallest k is
    returned with a warning.
    """
    if coph_df.empty:
        raise ValueError("Empty cophenetic dataframe — no valid NMF runs.")

    coph_df = coph_df.sort_values('k').reset_index(drop=True)

    # Drop k values where the NMF collapsed.
    # Two collapse signals:
    #   1. cophenetic = NaN  (C matrix was exactly degenerate → dist all-zero)
    #   2. mean_collapse_frac >= 0.90  (>90% of individuals in the dominant
    #      component across runs — consensus is stable but trivially so)
    valid = coph_df.dropna(subset=['cophenetic_correlation'])
    n_nan = len(coph_df) - len(valid)
    if n_nan > 0:
        nan_ks = coph_df.loc[coph_df['cophenetic_correlation'].isna(), 'k'].tolist()
        print(f"  WARNING: k={nan_ks} produced NaN cophenetic correlation "
              f"(NMF collapse — all individuals assigned to one component). "
              f"These k values are excluded from optimal-k selection.")

    if 'mean_collapse_frac' in valid.columns:
        collapsed_mask = valid['mean_collapse_frac'] >= 0.90
        if collapsed_mask.any():
            collapsed_ks = valid.loc[collapsed_mask, 'k'].tolist()
            print(f"  WARNING: k={collapsed_ks} excluded — NMF collapsed "
                  f"(≥90% individuals in dominant component despite non-NaN cophenetic). "
                  f"Consider increasing alpha_H or tightening feature filters.")
            valid = valid[~collapsed_mask]

    if valid.empty:
        fallback_k = int(coph_df['k'].min())
        print(f"  WARNING: All k values collapsed. Falling back to k={fallback_k}. "
              f"Results will not be meaningful — consider increasing alpha_H "
              f"(e.g. --alpha_H 0.5) or tightening feature filters.")
        return fallback_k

    coph_vals = valid['cophenetic_correlation'].values

    if method == 'max':
        best_idx = int(np.argmax(coph_vals))
    else:  # elbow
        if len(coph_vals) < 3:
            best_idx = int(np.argmax(coph_vals))
        else:
            # Second derivative (change in slope)
            d2 = np.diff(coph_vals, n=2)
            best_idx = int(np.argmin(d2)) + 1  # offset for double diff

    optimal_k = int(valid.iloc[best_idx]['k'])
    print(f"  Optimal k selected: {optimal_k} (method='{method}', "
          f"from {len(valid)} non-collapsed k values)")
    return optimal_k


def get_consensus_cluster_assignments(k_results, optimal_k, feature_names, index):
    """
    From the list of W matrices for the optimal k, build:
      - Average W matrix (soft membership)
      - Hard cluster assignments (argmax of averaged W rows)
      - Average H matrix (feature loadings)

    Returns
    -------
    assignments_df : pd.DataFrame  (IID index) with soft + hard assignments
    H_df           : pd.DataFrame  feature loadings per cluster
    """
    W_list = k_results[optimal_k]['W_list']
    H_list = k_results[optimal_k]['H_list']

    # Align components across runs (greedy column matching by correlation)
    W_ref = W_list[0]
    H_ref = H_list[0]
    W_aligned = [W_ref]
    H_aligned = [H_ref]

    for W, H in zip(W_list[1:], H_list[1:]):
        # Match columns of W to columns of W_ref using correlation
        corr_matrix = np.corrcoef(W.T, W_ref.T)[:optimal_k, optimal_k:]
        perm = []
        remaining = list(range(optimal_k))
        for i in range(optimal_k):
            best_j = max(remaining, key=lambda j: corr_matrix[i, j])
            perm.append(best_j)
            remaining.remove(best_j)
        W_aligned.append(W[:, perm])
        H_aligned.append(H[perm, :])

    W_avg = np.mean(W_aligned, axis=0)
    H_avg = np.mean(H_aligned, axis=0)

    # Normalize W rows to sum to 1 (soft membership fractions)
    row_sums = W_avg.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    W_norm = W_avg / row_sums

    cluster_labels = [f'cluster_{i+1}' for i in range(optimal_k)]

    assignments_df = pd.DataFrame(
        W_norm, index=index,
        columns=[f'membership_{c}' for c in cluster_labels]
    )
    assignments_df['hard_cluster'] = np.argmax(W_norm, axis=1) + 1
    assignments_df['cluster_label'] = [
        f'cluster_{i}' for i in assignments_df['hard_cluster']
    ]
    assignments_df['membership_entropy'] = -np.sum(
        np.where(W_norm > 0, W_norm * np.log(W_norm + 1e-12), 0), axis=1
    )

    H_df = pd.DataFrame(
        H_avg, index=cluster_labels, columns=feature_names
    )

    return assignments_df, H_df


# ============================================================================
# CLUSTER CHARACTERIZATION
# ============================================================================

def characterize_clusters(assignments_df, feature_matrix, phenotype_series=None):
    """
    Compute per-cluster summary statistics (median, IQR) for every feature
    in feature_matrix. Also includes phenotype prevalence if provided.

    Returns
    -------
    pd.DataFrame  rows = clusters, cols = features (median values)
    """
    merged = feature_matrix.join(
        assignments_df[['hard_cluster']], how='inner'
    )
    if phenotype_series is not None:
        merged = merged.join(phenotype_series.rename('PHENOTYPE'), how='left')

    summary_rows = []
    for cluster_id in sorted(merged['hard_cluster'].unique()):
        sub = merged[merged['hard_cluster'] == cluster_id]
        row = {'cluster': f'cluster_{cluster_id}', 'n': len(sub)}
        for col in feature_matrix.columns:
            if col in sub.columns:
                row[f'{col}_median'] = sub[col].median()
                row[f'{col}_iqr'] = sub[col].quantile(0.75) - sub[col].quantile(0.25)
        if 'PHENOTYPE' in sub.columns:
            n_cases = (sub['PHENOTYPE'] == 2).sum()
            row['case_prevalence'] = n_cases / len(sub) if len(sub) > 0 else np.nan
        summary_rows.append(row)

    return pd.DataFrame(summary_rows).set_index('cluster')


def map_clusters_to_concepts(H_df, feature_names):
    """
    Map feature loadings in H matrix to standardized clinical concepts
    using keyword matching for interpretability.

    Returns H_df with an additional summary of top features per cluster.
    """
    concept_map = {}
    for feat in feature_names:
        concept = map_feature_to_concept(str(feat))
        if concept:
            concept_map[feat] = concept

    top_features = {}
    for cluster in H_df.index:
        row = H_df.loc[cluster].sort_values(ascending=False)
        top_10 = row.head(10)
        top_features[cluster] = {
            'top_features': top_10.index.tolist(),
            'top_loadings': top_10.values.tolist(),
            'top_concepts': [concept_map.get(f, f) for f in top_10.index],
        }
    return top_features


# ============================================================================
# UDLER T2D VALIDATION
# ============================================================================

def compare_clusters_to_udler(cluster_profiles_df, clinical_cols):
    """
    Compare cluster clinical profiles to Udler T2D subtype reference profiles
    by computing cosine similarity of standardized median vectors.

    Parameters
    ----------
    cluster_profiles_df  rows=clusters, cols include {feature}_median
    clinical_cols        list of raw clinical feature names (before _median suffix)

    Returns
    -------
    similarity_df  pd.DataFrame  rows=clusters, cols=Udler subtypes, values=cosine sim
    best_match     dict          {cluster: best_matching_subtype}
    """
    # Map each clinical column to a concept
    col_to_concept = {}
    for col in clinical_cols:
        c = map_feature_to_concept(str(col))
        if c:
            col_to_concept[col] = c

    if not col_to_concept:
        print("  WARNING: no clinical features could be mapped to Udler concepts.")
        return pd.DataFrame(), {}

    # Build cluster vectors (z-score of medians across clusters)
    profile_data = {}
    for col, concept in col_to_concept.items():
        median_col = f'{col}_median'
        if median_col in cluster_profiles_df.columns:
            vals = cluster_profiles_df[median_col].values.astype(float)
            std = np.std(vals)
            if std > 0:
                profile_data[concept] = (vals - np.mean(vals)) / std
            else:
                profile_data[concept] = vals - np.mean(vals)

    if not profile_data:
        return pd.DataFrame(), {}

    profile_df = pd.DataFrame(profile_data, index=cluster_profiles_df.index)

    # Build Udler reference vectors
    all_concepts = profile_df.columns.tolist()
    similarity_rows = []
    for cluster in profile_df.index:
        cluster_vec = profile_df.loc[cluster].values
        row = {'cluster': cluster}
        for subtype, profile in UDLER_PROFILES.items():
            ref_vec = np.array([profile.get(c, 0) for c in all_concepts], dtype=float)
            norm_c = np.linalg.norm(cluster_vec)
            norm_r = np.linalg.norm(ref_vec)
            if norm_c > 0 and norm_r > 0:
                sim = float(np.dot(cluster_vec, ref_vec) / (norm_c * norm_r))
            else:
                sim = 0.0
            row[subtype] = sim
        similarity_rows.append(row)

    sim_df = pd.DataFrame(similarity_rows).set_index('cluster')
    best_match = sim_df.idxmax(axis=1).to_dict()

    print("\n  Cluster → Udler subtype similarity:")
    for cluster, match in best_match.items():
        print(f"    {cluster} → {match} "
              f"(sim={sim_df.loc[cluster, match]:.3f})")

    return sim_df, best_match


# ============================================================================
# VISUALIZATIONS
# ============================================================================

def plot_cophenetic_curve(coph_df, output_path):
    """Plot cophenetic correlation and reconstruction error vs k."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    ax.plot(coph_df['k'], coph_df['cophenetic_correlation'],
            'o-', color='#0072B2', linewidth=2, markersize=7)
    ax.set_xlabel('Number of clusters (k)', fontsize=11)
    ax.set_ylabel('Cophenetic correlation', fontsize=11)
    ax.set_title('Cluster stability (cophenetic correlation)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(coph_df['k'])

    ax2 = axes[1]
    ax2.plot(coph_df['k'], coph_df['mean_reconstruction_error'],
             's--', color='#D55E00', linewidth=2, markersize=7)
    ax2.set_xlabel('Number of clusters (k)', fontsize=11)
    ax2.set_ylabel('Mean reconstruction error', fontsize=11)
    ax2.set_title('NMF reconstruction error', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(coph_df['k'])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved cophenetic curve → {output_path}")


def plot_consensus_map(C, k, output_path):
    """Heatmap of consensus matrix for a given k."""
    fig, ax = plt.subplots(figsize=(8, 7))
    # Sort by hierarchical clustering for visual clarity
    dist = 1 - C
    np.fill_diagonal(dist, 0)
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method='average')
    order = fcluster(Z, k, criterion='maxclust')
    idx = np.argsort(order)
    C_sorted = C[np.ix_(idx, idx)]
    sns.heatmap(
        C_sorted, ax=ax,
        cmap='Blues', vmin=0, vmax=1,
        xticklabels=False, yticklabels=False,
        cbar_kws={'label': 'Consensus co-clustering probability'},
    )
    ax.set_title(f'Consensus matrix (k={k})', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved consensus map (k={k}) → {output_path}")


def plot_cluster_profiles(H_df, output_path, top_n=20):
    """
    Heatmap of top features per cluster from H matrix (normalized loadings).
    """
    # Select top_n features by max loading across clusters
    max_loading = H_df.max(axis=0)
    top_features = max_loading.nlargest(top_n).index.tolist()
    H_sub = H_df[top_features]

    # Column-normalize for relative comparison
    col_max = H_sub.max(axis=0)
    col_max[col_max == 0] = 1
    H_norm = H_sub / col_max

    fig, ax = plt.subplots(figsize=(max(8, top_n * 0.4), max(4, len(H_df) * 0.8)))
    sns.heatmap(
        H_norm, ax=ax,
        cmap='YlOrRd', vmin=0, vmax=1,
        linewidths=0.3, linecolor='#dddddd',
        cbar_kws={'label': 'Normalized feature loading'},
    )
    ax.set_title('Top feature loadings per cluster (H matrix)', fontsize=12)
    ax.set_xlabel('Features', fontsize=10)
    ax.set_ylabel('Cluster', fontsize=10)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved cluster profiles → {output_path}")


def plot_cluster_clinical_boxplots(assignments_df, clinical_matrix,
                                   clinical_cols_of_interest, output_path,
                                   udler_match=None):
    """
    Boxplots of key clinical features stratified by cluster assignment.
    """
    merged = clinical_matrix.join(assignments_df[['cluster_label']], how='inner')
    if merged.empty:
        return

    # Filter to available columns
    plot_cols = [c for c in clinical_cols_of_interest if c in merged.columns]
    if not plot_cols:
        plot_cols = [c for c in merged.columns if c != 'cluster_label'][:8]

    n_cols = min(4, len(plot_cols))
    n_rows = int(np.ceil(len(plot_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 4, n_rows * 3.5))
    axes = np.array(axes).flatten()

    clusters = sorted(merged['cluster_label'].unique())
    palette = {c: CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
               for i, c in enumerate(clusters)}

    for ax, col in zip(axes, plot_cols):
        data = [merged.loc[merged['cluster_label'] == c, col].dropna()
                for c in clusters]
        bp = ax.boxplot(data, patch_artist=True, notch=False,
                        medianprops={'color': 'black', 'linewidth': 2})
        for patch, c in zip(bp['boxes'], clusters):
            patch.set_facecolor(palette[c])
            patch.set_alpha(0.75)

        ax.set_xticklabels(clusters, rotation=30, ha='right', fontsize=8)
        label = col
        if udler_match:
            concept = map_feature_to_concept(col)
        ax.set_title(label, fontsize=9)
        ax.set_ylabel('Value', fontsize=8)
        ax.grid(True, axis='y', alpha=0.3)

    # Hide unused axes
    for ax in axes[len(plot_cols):]:
        ax.set_visible(False)

    # Legend for Udler match
    if udler_match:
        legend_patches = [
            mpatches.Patch(color=palette.get(c, '#999999'),
                           label=f"{c} → {udler_match.get(c, '?')}")
            for c in clusters
        ]
        fig.legend(handles=legend_patches, loc='lower right',
                   title='Cluster (→ Udler subtype)', fontsize=8,
                   title_fontsize=9)

    plt.suptitle('Clinical feature distributions by cluster', fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved clinical boxplots → {output_path}")


def plot_udler_similarity(sim_df, best_match, output_path):
    """Heatmap of cluster-to-Udler subtype cosine similarities."""
    if sim_df.empty:
        return
    fig, ax = plt.subplots(figsize=(7, max(3, len(sim_df) * 0.8)))
    sns.heatmap(
        sim_df, ax=ax, annot=True, fmt='.2f',
        cmap='RdYlGn', vmin=-1, vmax=1, center=0,
        linewidths=0.5, linecolor='#cccccc',
        cbar_kws={'label': 'Cosine similarity'},
    )
    ax.set_title('Cluster similarity to Udler T2D subtypes', fontsize=12)
    ax.set_xlabel('Udler subtype', fontsize=10)
    ax.set_ylabel('Cluster', fontsize=10)

    # Annotate best match
    for i, cluster in enumerate(sim_df.index):
        best = best_match.get(cluster, '')
        if best in sim_df.columns:
            j = list(sim_df.columns).index(best)
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                       edgecolor='black', lw=2))
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved Udler similarity plot → {output_path}")


def plot_prs_cluster_distributions(assignments_df, prs_matrix, output_path):
    """Violin plots of PRS model scores stratified by cluster."""
    merged = prs_matrix.join(assignments_df[['cluster_label']], how='inner')
    if merged.empty:
        return

    prs_cols = [c for c in prs_matrix.columns]
    n_cols = min(3, len(prs_cols))
    n_rows = int(np.ceil(len(prs_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 4.5, n_rows * 3.5))
    axes = np.array(axes).flatten()

    clusters = sorted(merged['cluster_label'].unique())
    palette = {c: CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
               for i, c in enumerate(clusters)}

    for ax, col in zip(axes, prs_cols):
        plot_data = merged[[col, 'cluster_label']].dropna()
        if plot_data.empty:
            ax.set_visible(False)
            continue
        sns.violinplot(
            data=plot_data, x='cluster_label', y=col,
            palette=palette, ax=ax, inner='quartile',
            order=clusters,
        )
        ax.set_title(col.replace('_prs', '').replace('_', ' ').title(),
                     fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('Scaled PRS', fontsize=9)
        ax.tick_params(axis='x', rotation=30)
        ax.grid(True, axis='y', alpha=0.3)

    for ax in axes[len(prs_cols):]:
        ax.set_visible(False)

    plt.suptitle('PRS model score distributions by cluster', fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved PRS cluster distributions → {output_path}")


# ============================================================================
# MAIN
# ============================================================================

def main(
    pheno_data,
    clinical_file,
    mode,
    k_min,
    k_max,
    n_runs,
    l1_ratio,
    alpha_w,
    k_select_method,
    phenotype_name,
    important_features_file,
    cases_only,
    max_iter,
):
    scores_path = os.path.join(pheno_data, 'scores')
    out_dir = os.path.join(scores_path, 'bnmf')
    fig_dir = os.path.join(pheno_data, 'figures', 'bnmf')
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  prs_bnmf_subtype_analysis")
    print(f"  phenotype : {phenotype_name}")
    print(f"  mode      : {mode}")
    print(f"  k range   : {k_min} – {k_max}  ({n_runs} runs each)")
    print(f"{'='*60}\n")

    # ------------------------------------------------------------------
    # 1. Load PRS files
    # ------------------------------------------------------------------
    print("[1] Loading PRS files ...")
    # Detect available models from combined PRS file if present
    combined_path = os.path.join(scores_path, 'combinedPRSGroups.filtered.csv')
    models_to_use = None
    if os.path.exists(combined_path):
        combined_df = pd.read_csv(combined_path)
        bin_cols = [c for c in combined_df.columns if c.startswith('bin_')]
        models_to_use = [c.replace('bin_', '') for c in bin_cols]
        print(f"  Models from combined file: {models_to_use}")

    prs_dict = load_prs_files(scores_path, models_to_use)
    if not prs_dict:
        raise FileNotFoundError(
            f"No *.mixed.prs.csv files found under {scores_path}/prsScores/. "
            "Run calculate_prs_post_modelling.py first."
        )

    # ------------------------------------------------------------------
    # 2. Build feature matrix
    # ------------------------------------------------------------------
    print(f"\n[2] Building feature matrix (mode='{mode}') ...")
    important_feats_df = None
    if important_features_file and os.path.exists(important_features_file):
        important_feats_df = pd.read_csv(important_features_file)
        print(f"  Loaded important features: {len(important_feats_df)} rows")

    if mode == 'prs_scores':
        prs_matrix, phenotype_series = build_prs_score_matrix(prs_dict)
    else:  # weighted_features
        prs_matrix, phenotype_series = build_weighted_feature_matrix(
            prs_dict, important_feats_df
        )

    print(f"  PRS/feature matrix: {prs_matrix.shape}")

    # ------------------------------------------------------------------
    # 3. Load clinical measures
    # ------------------------------------------------------------------
    print("\n[3] Loading clinical measures ...")
    clinical_matrix = load_clinical_measures(clinical_file)

    # ------------------------------------------------------------------
    # 4. Filter individuals
    # ------------------------------------------------------------------
    shared_iids = prs_matrix.index.intersection(clinical_matrix.index)
    print(f"\n[4] Shared individuals (PRS ∩ clinical): {len(shared_iids)}")

    if cases_only and phenotype_series is not None:
        case_iids = phenotype_series[phenotype_series == 2].index
        shared_iids = shared_iids.intersection(case_iids)
        print(f"  Restricted to cases (PHENOTYPE==2): {len(shared_iids)}")

    if len(shared_iids) < 50:
        raise ValueError(
            f"Too few individuals for clustering ({len(shared_iids)}). "
            "Check that IIDs match between PRS and clinical files."
        )

    prs_sub = prs_matrix.loc[shared_iids]
    clin_sub = clinical_matrix.loc[shared_iids]

    # ------------------------------------------------------------------
    # 5. Association analysis: PRS models ↔ clinical measures
    # ------------------------------------------------------------------
    print("\n[5] PRS–clinical association analysis ...")
    if mode == 'prs_scores':
        rho_pivot, assoc_df = compute_prs_clinical_associations(prs_sub, clin_sub)
    else:
        # For weighted_features mode, use a summarized version per model
        # (mean of feature weights per model) to keep the association analysis tractable
        model_means = {}
        for model_name, df in prs_dict.items():
            if 'scaled_prs' in df.columns:
                model_means[f'{model_name}_prs'] = df['scaled_prs']
            elif 'prs' in df.columns:
                model_means[f'{model_name}_prs'] = df['prs']
        if model_means:
            prs_scores_for_assoc = pd.DataFrame(model_means).loc[
                prs_matrix.index.intersection(list(model_means.values())[0].index)
            ]
            prs_scores_for_assoc = prs_scores_for_assoc.loc[shared_iids]
            rho_pivot, assoc_df = compute_prs_clinical_associations(
                prs_scores_for_assoc, clin_sub
            )
        else:
            rho_pivot, assoc_df = pd.DataFrame(), pd.DataFrame()

    if not assoc_df.empty:
        assoc_df.to_csv(
            os.path.join(out_dir, 'prs_clinical_associations.csv'), index=False
        )
        plot_association_heatmap(
            rho_pivot,
            os.path.join(fig_dir, 'prs_clinical_heatmap.png'),
            title=f'{phenotype_name}: PRS–Clinical Spearman ρ'
        )

    # ------------------------------------------------------------------
    # 6. Prepare NMF input matrix
    # ------------------------------------------------------------------
    print("\n[6] Preparing NMF input matrix ...")

    if mode == 'prs_scores':
        # Combine PRS scores + clinical measures
        combined = prs_sub.join(clin_sub, rsuffix='_clin', how='inner')
    else:
        # Combine weighted features + clinical measures
        combined = prs_sub.loc[shared_iids].join(
            clin_sub, rsuffix='_clin', how='inner'
        )

    print(f"  Combined matrix before preprocessing: {combined.shape}")
    X_scaled, retained_cols = impute_and_scale(combined)
    print(f"  After imputation/scaling: {X_scaled.shape}")

    # Cast to float32 to halve NMF working memory (sklearn NMF accepts float32).
    X_arr = X_scaled.values.astype(np.float32)
    n_samples, n_features = X_arr.shape
    print(f"  NMF input: {n_samples} samples × {n_features} features  "
          f"({X_arr.nbytes / 1e6:.1f} MB, float32)")

    # Add small epsilon to avoid zero rows/cols (NMF requirement)
    X_arr = X_arr + np.float32(1e-9)

    # ------------------------------------------------------------------
    # 7. Consensus bNMF
    # ------------------------------------------------------------------
    print(f"\n[7] Running consensus bNMF (k={k_min}–{k_max}, {n_runs} runs) ...")
    k_results, coph_df = run_consensus_bnmf(
        X_arr, k_min=k_min, k_max=k_max,
        n_runs=n_runs, l1_ratio=l1_ratio, alpha_W=alpha_w,
        max_iter=max_iter, verbose=True,
    )

    coph_df.to_csv(os.path.join(out_dir, 'cophenetic_curve.csv'), index=False)
    plot_cophenetic_curve(coph_df, os.path.join(fig_dir, 'cophenetic_curve.png'))

    # ------------------------------------------------------------------
    # 8. Select optimal k  (done here, before plotting, so we know which
    #    k's W_list must be kept and which can be freed on the fly)
    # ------------------------------------------------------------------
    print(f"\n[8] Selecting optimal k ({k_select_method} method) ...")
    optimal_k = select_optimal_k(coph_df, method=k_select_method)

    # Plot consensus maps by recomputing C one k at a time so we never hold
    # more than one n×n matrix in memory at once.
    # For large cohorts subsample for visualisation — the full n×n matrix is
    # unreadable at n>10k and would be tens of GB (67k × 67k × 4 B ≈ 18 GB).
    _rng = np.random.default_rng(42)
    for k in sorted(k_results.keys()):
        W_list_k = k_results[k]['W_list']
        n_full = W_list_k[0].shape[0]
        if n_full > MAX_CONSENSUS_PLOT_N:
            idx = _rng.choice(n_full, MAX_CONSENSUS_PLOT_N, replace=False)
            idx.sort()
            W_list_plot = [W[idx] for W in W_list_k]
            print(f"  Consensus map k={k}: subsampled {MAX_CONSENSUS_PLOT_N} "
                  f"of {n_full} individuals for visualisation")
        else:
            W_list_plot = W_list_k
        C_plot = _soft_consensus_matrix(W_list_plot)
        plot_consensus_map(
            C_plot, k,
            os.path.join(fig_dir, f'consensus_map_k{k}.png')
        )
        del C_plot, W_list_plot
        # Free W/H for non-optimal k immediately — no longer needed.
        if k != optimal_k:
            del k_results[k]['W_list'], k_results[k]['H_list']
        gc.collect()

    assignments_df, H_df = get_consensus_cluster_assignments(
        k_results, optimal_k,
        feature_names=retained_cols,
        index=X_scaled.index,
    )

    # Add phenotype column if available
    if phenotype_series is not None:
        assignments_df = assignments_df.join(
            phenotype_series.rename('PHENOTYPE'), how='left'
        )

    assignments_df.to_csv(
        os.path.join(out_dir, 'cluster_assignments.csv')
    )
    H_df.to_csv(os.path.join(out_dir, 'feature_loadings.csv'))
    print(f"  Saved cluster assignments ({optimal_k} clusters) → {out_dir}")

    # Cluster size summary
    size_summary = assignments_df['hard_cluster'].value_counts().sort_index()
    print("\n  Cluster sizes:")
    for cl, n in size_summary.items():
        print(f"    cluster_{cl}: n={n}")

    # ------------------------------------------------------------------
    # 9. Cluster profiles using clinical measures
    # ------------------------------------------------------------------
    print("\n[9] Characterizing clusters with clinical measures ...")
    cluster_profiles = characterize_clusters(
        assignments_df, clin_sub, phenotype_series
    )
    cluster_profiles.to_csv(os.path.join(out_dir, 'cluster_clinical_profiles.csv'))
    print(f"  Saved cluster clinical profiles → {out_dir}")

    # Feature loadings profile
    plot_cluster_profiles(
        H_df, os.path.join(fig_dir, 'cluster_feature_loadings.png'), top_n=30
    )

    # Clinical boxplots (use columns with identifiable clinical keywords)
    clinical_interest = [
        c for c in clin_sub.columns
        if map_feature_to_concept(str(c)) is not None
    ]
    if not clinical_interest:
        clinical_interest = list(clin_sub.columns[:12])

    # ------------------------------------------------------------------
    # 10. T2D Udler validation (if applicable)
    # ------------------------------------------------------------------
    udler_match = None
    if 'diabetes' in phenotype_name.lower() or 't2d' in phenotype_name.lower():
        print("\n[10] Comparing clusters to Udler T2D subtypes ...")
        sim_df, udler_match = compare_clusters_to_udler(
            cluster_profiles,
            [c for c in clin_sub.columns]
        )
        if not sim_df.empty:
            sim_df.to_csv(os.path.join(out_dir, 'udler_similarity.csv'))
            plot_udler_similarity(
                sim_df, udler_match,
                os.path.join(fig_dir, 'cluster_udler_comparison.png')
            )
            # Add Udler match to assignments
            assignments_df['udler_match'] = assignments_df['cluster_label'].map(
                udler_match
            )
            assignments_df.to_csv(
                os.path.join(out_dir, 'cluster_assignments.csv')
            )

    plot_cluster_clinical_boxplots(
        assignments_df, clin_sub,
        clinical_interest,
        os.path.join(fig_dir, 'cluster_clinical_boxplots.png'),
        udler_match=udler_match,
    )

    # PRS distributions per cluster (only for prs_scores mode)
    if mode == 'prs_scores':
        plot_prs_cluster_distributions(
            assignments_df, prs_sub,
            os.path.join(fig_dir, 'prs_cluster_distributions.png'),
        )

    # ------------------------------------------------------------------
    # 11. Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  ANALYSIS COMPLETE")
    print(f"  Optimal k : {optimal_k}")
    print(f"  Outputs   : {out_dir}")
    print(f"  Figures   : {fig_dir}")
    print(f"{'='*60}\n")

    # Top feature interpretation per cluster
    top_feats = map_clusters_to_concepts(H_df, retained_cols)
    print("  Top features per cluster:")
    for cluster, info in top_feats.items():
        print(f"    {cluster}: {info['top_concepts'][:5]}")

    return assignments_df, H_df, coph_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='bNMF subtype analysis of ePRS models and clinical measures.',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
      '--pheno_data', required=True,
      help='Path to phenotype results directory (e.g. results/type2Diabetes/summedEpi)'
    )
    parser.add_argument(
      '--clinical_file', required=True,
      help='Path to clinical measures CSV (participant_environment.csv or similar)'
    )
    parser.add_argument(
        '--mode', default='prs_scores',
        choices=['prs_scores', 'weighted_features'],
        help=(
            'Feature input mode:\n'
            '  prs_scores        – aggregate scaled_prs per model + clinical measures\n'
            '  weighted_features – per-individual weighted feature contributions + clinical measures'
        )
    )
    parser.add_argument(
        '--k_min', type=int, default=2,
        help='Minimum number of clusters to test (default: 2)'
    )
    parser.add_argument(
        '--k_max', type=int, default=8,
        help='Maximum number of clusters to test (default: 8)'
    )
    parser.add_argument(
        '--n_runs', type=int, default=30,
        help='Number of NMF restarts per k for consensus (default: 30)'
    )
    parser.add_argument(
        '--l1_ratio', type=float, default=0.1,
        help='L1 regularization ratio for ARD-like sparsity on W (default: 0.1)'
    )
    parser.add_argument(
        '--alpha_w', type=float, default=0.1,
        help='Regularization strength on W (default: 0.1)'
    )
    parser.add_argument(
        '--k_select_method', default='elbow',
        choices=['elbow', 'max'],
        help='Method for selecting optimal k from cophenetic curve (default: elbow)'
    )
    parser.add_argument(
        '--phenotype_name', default='phenotype',
        help='Human-readable phenotype label used in plot titles (default: phenotype)'
    )
    parser.add_argument(
        '--important_features_file', default=None,
        help=(
            'Path to important features CSV (e.g. importantFeaturesPostShap.csv).\n'
            'Used in weighted_features mode to restrict to SHAP-selected features.'
        )
    )
    parser.add_argument(
        '--cases_only', action='store_true',
        help='Restrict bNMF to cases (PHENOTYPE==2) only'
    )
    parser.add_argument(
        '--max_iter', type=int, default=500,
        help='Maximum NMF iterations per run (default: 500)'
    )

    args = parser.parse_args()
    
#   phenotype_name = 'type2Diabetes'
#   pheno_data = f'/Users/kerimulterer/prsInteractive/results/{phenotype_name}/combinedAnalysis'
#   clinical_file = '/Users/kerimulterer/prsInteractive/results/participant_environment.csv'
#   important_features_file = f'/Users/kerimulterer/prsInteractive/results/{phenotype_name}/combinedAnalysis/scores/featureScoresReducedFinalModel.filtered.csv'

    pheno_data = args.pheno_data or os.environ.get('PHENO_DATA')
    if not pheno_data:
      raise ValueError("Provide --pheno_data or set PHENO_DATA env var.")

    main(
        pheno_data=pheno_data,
        clinical_file=args.clinical_file,
        mode=args.mode,
        k_min=args.k_min,
        k_max=args.k_max,
        n_runs=args.n_runs,
        l1_ratio=args.l1_ratio,
        alpha_w=args.alpha_w,
        k_select_method=args.k_select_method,
        phenotype_name=args.phenotype_name,
        important_features_file=args.important_features_file,
        cases_only=args.cases_only,
        max_iter=args.max_iter,
    )
  
# main(
#   pheno_data=pheno_data,
#   clinical_file=args.clinical_file,
#   mode=args.mode,
#   k_min=args.k_min,
#   k_max=args.k_max,
#   n_runs=args.n_runs,
#   l1_ratio=args.l1_ratio,
#   alpha_w=args.alpha_w,
#   k_select_method=args.k_select_method,
#   phenotype_name=args.phenotype_name,
#   important_features_file=args.important_features_file,
#   cases_only=args.cases_only,
#   max_iter=args.max_iter,
# )