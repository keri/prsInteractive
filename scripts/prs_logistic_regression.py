#!/usr/bin/env python3

# Keri Multerer 2025
# Unpenalized logistic regression:  PRSComp + covariates  vs
#                                   PRSComp + covariates + clinical measure
# 5-fold stratified cross-validation on the holdout set.
#
# Scaling note
# ------------
# The PRS (scaled_prs_*) columns are already standardised upstream.
# Clinical measures (HbA1c, BMI, Glucose) and continuous covariates
# (age, PC1-PC10) span very different scales, so they are z-scored per
# training fold to avoid numerical instability in the unpenalized solver.
# SEX is binary and left unscaled.

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, roc_curve
)
from scipy.stats import norm

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _count_reclassifications(old_risk, new_risk, is_case):
    """Pencina 2008 reclassification counter (mirrors clinical_measure_performance)."""
    is_control = ~is_case
    up_cases      = int(np.sum(~old_risk[is_case]    &  new_risk[is_case]))
    down_cases    = int(np.sum( old_risk[is_case]    & ~new_risk[is_case]))
    up_controls   = int(np.sum(~old_risk[is_control] &  new_risk[is_control]))
    down_controls = int(np.sum( old_risk[is_control] & ~new_risk[is_control]))
    return {
        'up_cases': up_cases, 'down_cases': down_cases,
        'up_controls': up_controls, 'down_controls': down_controls,
        'n_cases': int(is_case.sum()), 'n_controls': int(is_control.sum()),
    }


def calculate_nri_from_arrays(y_true, pred_ref, pred_new, threshold=0.5):
    """
    NRI with 95% CI comparing two sets of predicted probabilities.

    pred_ref / pred_new are binarised at `threshold` (default 0.5).
    For the holdout CV context we use the Youden-optimal threshold derived
    from the training folds; pass it in via `threshold`.

    Returns
    -------
    dict with keys: nri, nri_events, nri_non_events, ci_low, ci_high, se
    """
    old_risk = pred_ref >= threshold
    new_risk = pred_new >= threshold
    is_case  = y_true.astype(bool)

    counts = _count_reclassifications(old_risk, new_risk, is_case)
    n_cases    = counts['n_cases']
    n_controls = counts['n_controls']

    nri_events     = (counts['up_cases']      - counts['down_cases'])    / n_cases
    nri_non_events = (counts['down_controls'] - counts['up_controls'])   / n_controls
    nri            = nri_events + nri_non_events

    se_events     = np.sqrt(counts['up_cases']    + counts['down_cases'])    / n_cases
    se_non_events = np.sqrt(counts['up_controls'] + counts['down_controls']) / n_controls
    se_nri        = np.sqrt(se_events**2 + se_non_events**2)

    ci_low  = nri - 1.96 * se_nri
    ci_high = nri + 1.96 * se_nri

    return {
        'nri':            nri,
        'nri_events':     nri_events,
        'nri_non_events': nri_non_events,
        'ci_low':         ci_low,
        'ci_high':        ci_high,
        'se':             se_nri,
        'n_cases':        n_cases,
        'n_controls':     n_controls,
    }


def youden_threshold(y_true, y_prob):
    """Return the probability threshold maximising sensitivity + specificity."""
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    idx = np.argmax(tpr - fpr)
    return float(thresholds[idx])


# ---------------------------------------------------------------------------
# core analysis
# ---------------------------------------------------------------------------

def run_logistic_cv(df, prs_col, covariate_cols, clinical_cols,
                    outcome_col='PHENOTYPE', n_splits=5, seed=42):
    """
    5-fold stratified CV comparing:
      Model 1 (base)    : PRSComp + covariates
      Model 2 (augmented): PRSComp + covariates + one clinical measure

    Features that need scaling (continuous): covariates + clinical measures.
    SEX is left unscaled; PRS is already scaled.

    Returns
    -------
    dict keyed by clinical measure name, each containing per-fold and
    pooled metrics: AUC, precision, recall, NRI (with 95% CI).
    Also returns base model pooled metrics under key '_base'.
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

    # Columns that should be z-scored on each fold
    cols_to_scale = [c for c in covariate_cols if c != 'SEX']

    # Full feature matrix (base model, unscaled continuous)
    base_features = [prs_col] + covariate_cols
    y = df[outcome_col].values

    # Storage: fold-level predictions for base model
    base_prob_pool = np.full(len(df), np.nan)
    base_thresh_per_fold = []

    # fold-level predictions per clinical measure
    clin_prob_pool = {c: np.full(len(df), np.nan) for c in clinical_cols}

    indices = np.arange(len(df))

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(indices, y)):
        df_train = df.iloc[train_idx].copy()
        df_test  = df.iloc[test_idx].copy()

        # Fit scaler on training fold continuous features (excluding PRS & SEX)
        scaler = StandardScaler()
        scaler.fit(df_train[cols_to_scale])

        df_train[cols_to_scale] = scaler.transform(df_train[cols_to_scale])
        df_test[cols_to_scale]  = scaler.transform(df_test[cols_to_scale])

        y_train = df_train[outcome_col].values
        y_test  = df_test[outcome_col].values

        # --- Base model ---
        X_train_base = df_train[base_features].values
        X_test_base  = df_test[base_features].values

        base_model = LogisticRegression(penalty=None, solver='lbfgs',
                                        max_iter=2000, random_state=seed)
        base_model.fit(X_train_base, y_train)
        base_prob_fold = base_model.predict_proba(X_test_base)[:, 1]
        base_prob_pool[test_idx] = base_prob_fold

        # Youden threshold from training fold predictions
        train_prob_base = base_model.predict_proba(X_train_base)[:, 1]
        thresh = youden_threshold(y_train, train_prob_base)
        base_thresh_per_fold.append(thresh)

        # --- Augmented model: one clinical measure at a time ---
        for clin_col in clinical_cols:
            # Scale the clinical feature on this fold
            clin_scaler = StandardScaler()
            clin_train = clin_scaler.fit_transform(df_train[[clin_col]])
            clin_test  = clin_scaler.transform(df_test[[clin_col]])

            X_train_aug = np.hstack([X_train_base, clin_train])
            X_test_aug  = np.hstack([X_test_base,  clin_test])

            aug_model = LogisticRegression(penalty=None, solver='lbfgs',
                                           max_iter=2000, random_state=seed)
            try:
                aug_model.fit(X_train_aug, y_train)
                aug_prob_fold = aug_model.predict_proba(X_test_aug)[:, 1]
            except Exception as exc:
                print(f'  [WARN] fold {fold_idx} {clin_col}: {exc}')
                aug_prob_fold = base_prob_fold  # fall back

            clin_prob_pool[clin_col][test_idx] = aug_prob_fold

    # -------------------------------------------------------------------
    # Pool predictions and compute metrics
    # -------------------------------------------------------------------
    # Use median Youden threshold from training folds for NRI binarisation
    base_thresh = float(np.median(base_thresh_per_fold))
    print(f'\nMedian Youden threshold (base model): {base_thresh:.3f}')

    def pooled_metrics(y_true, prob_ref, prob_new, thresh_ref):
        auc    = roc_auc_score(y_true, prob_new)
        preds  = (prob_new >= thresh_ref).astype(int)
        prec   = precision_score(y_true, preds, zero_division=0)
        rec    = recall_score(y_true, preds, zero_division=0)
        nri    = calculate_nri_from_arrays(y_true, prob_ref, prob_new,
                                           threshold=thresh_ref)
        return {'auc': auc, 'precision': prec, 'recall': rec, **nri}

    results = {}

    # Base model metrics (NRI self-comparison → should be 0, used for reference AUC)
    base_auc       = roc_auc_score(y, base_prob_pool)
    base_preds     = (base_prob_pool >= base_thresh).astype(int)
    base_precision = precision_score(y, base_preds, zero_division=0)
    base_recall    = recall_score(y, base_preds, zero_division=0)
    results['_base'] = {
        'auc':       base_auc,
        'precision': base_precision,
        'recall':    base_recall,
    }
    print(f'\n--- Base model (PRSComp + covariates) ---')
    print(f'  AUC={base_auc:.4f}  Precision={base_precision:.4f}  Recall={base_recall:.4f}')

    for clin_col in clinical_cols:
        m = pooled_metrics(y, base_prob_pool, clin_prob_pool[clin_col], base_thresh)
        results[clin_col] = m
        print(f'\n--- Augmented: + {clin_col} ---')
        print(f'  AUC={m["auc"]:.4f}  Precision={m["precision"]:.4f}  Recall={m["recall"]:.4f}')
        print(f'  NRI={m["nri"]:.4f}  95%CI=[{m["ci_low"]:.4f}, {m["ci_high"]:.4f}]')
        print(f'  NRI events={m["nri_events"]:.4f}  NRI non-events={m["nri_non_events"]:.4f}')

    return results, base_prob_pool, clin_prob_pool


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_auc_comparison(results, base_auc, figPath, prs_col_label='PRSComp'):
    """
    Horizontal bar chart: AUC for PRSComp+covariate vs PRSComp alone.

    Bars represent the augmented model AUC; the base model AUC is shown
    as a vertical dashed reference line.
    """
    clin_keys  = [k for k in results if k != '_base']
    labels     = [k.replace('_', ' ') for k in clin_keys]
    aucs       = [results[k]['auc'] for k in clin_keys]

    fig, ax = plt.subplots(figsize=(7, 0.7 * len(clin_keys) + 2))

    colors = ['#BA7517' if a > base_auc else '#378ADD' for a in aucs]

    y_pos = np.arange(len(labels))
    ax.barh(y_pos, aucs, color=colors, height=0.5, zorder=2)

    ax.axvline(base_auc, color='#2C2C2A', linewidth=1.4, linestyle='--',
               label=f'{prs_col_label} alone (AUC={base_auc:.3f})')

    for i, (a, label) in enumerate(zip(aucs, labels)):
        ax.text(a + 0.002, i, f'{a:.3f}', va='center', fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('AUC (5-fold CV)', fontsize=11)
    ax.set_title(f'AUC: {prs_col_label} + clinical measure vs {prs_col_label} alone',
                 fontsize=12, pad=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # x-axis range: start just below minimum AUC
    all_aucs = aucs + [base_auc]
    ax.set_xlim(max(0, min(all_aucs) - 0.05), min(1, max(all_aucs) + 0.06))

    fig.tight_layout()
    out = os.path.join(figPath, 'auc_prs_vs_prs_plus_clinical.png')
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'\nAUC comparison plot saved to {out}')
    return fig


def plot_nri_forest(results, figPath):
    """
    Forest-style plot of NRI with 95% CI for each clinical measure augmentation.
    """
    clin_keys = [k for k in results if k != '_base']
    labels    = [k.replace('_', ' ') for k in clin_keys]
    nris      = [results[k]['nri']     for k in clin_keys]
    ci_los    = [results[k]['ci_low']  for k in clin_keys]
    ci_his    = [results[k]['ci_high'] for k in clin_keys]

    fig, ax = plt.subplots(figsize=(7, 0.7 * len(clin_keys) + 2))

    y_pos = np.arange(len(labels))
    xerr_lo = np.array(nris) - np.array(ci_los)
    xerr_hi = np.array(ci_his) - np.array(nris)

    colors = ['#BA7517' if n > 0 else '#378ADD' for n in nris]
    ax.errorbar(nris, y_pos, xerr=[xerr_lo, xerr_hi],
                fmt='D', color='#2C2C2A', markersize=7, capsize=4,
                linewidth=1.2, zorder=3)
    for i, (n, c) in enumerate(zip(nris, colors)):
        ax.scatter(n, i, color=c, s=60, zorder=4)

    ax.axvline(0, color='#2C2C2A', linewidth=1, linestyle='--', alpha=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('NRI vs PRSComp alone (95% CI)', fontsize=11)
    ax.set_title('Net Reclassification Improvement\n(PRSComp + covariate + clinical measure vs PRSComp + covariate)',
                 fontsize=11, pad=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()

    out = os.path.join(figPath, 'nri_forest_prs_plus_clinical.png')
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'NRI forest plot saved to {out}')
    return fig


# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------

def load_data(pheno_data, results_path, pheno, prs_col,
              covariate_cols, clinical_source_path=None):
    """
    Load and merge holdout scores, participant covariates, and clinical data.

    Returns
    -------
    pd.DataFrame with columns: IID, PHENOTYPE, prs_col, covariate_cols,
                                + clinical columns
    list : clinical column names present in the merged frame
    """
    # PRS + PHENOTYPE (already scaled)
    holdout_path = os.path.join(pheno_data, 'scores', 'combinedPRSGroups.holdout.filtered.csv')
    df_prs = pd.read_csv(holdout_path)

    # Recode phenotype: UKB convention 1=control,2=case → 0/1
    df_prs['PHENOTYPE'] = df_prs['PHENOTYPE'].replace({1: 0, 2: 1})

    if prs_col not in df_prs.columns:
        available = [c for c in df_prs.columns if 'scaled_prs' in c]
        raise ValueError(
            f"PRS column '{prs_col}' not found. Available scaled_prs columns: {available}"
        )

    # Participant covariates (age, SEX, PCs)
    env_path = os.path.join(results_path, 'participant_environment.csv')
    keep_env = ['Participant ID'] + covariate_cols
    keep_env_present = [c for c in keep_env if c in pd.read_csv(env_path, nrows=0).columns]
    df_env = pd.read_csv(env_path, usecols=keep_env_present)
    df_env.rename(columns={'Participant ID': 'IID'}, inplace=True)

    missing_cov = [c for c in covariate_cols if c not in df_env.columns]
    if missing_cov:
        print(f'[WARN] Covariate columns not found in participant_environment.csv: {missing_cov}')

    # Clinical measures
    if clinical_source_path is None:
        parent = os.path.dirname(os.path.abspath(pheno_data))
        for cand in ('productEpi', 'summedEpi'):
            fp = os.path.join(parent, cand, 'clinicalEnvironmentHoldout.csv')
            if os.path.exists(fp):
                clinical_source_path = fp
                break
        if clinical_source_path is None:
            raise FileNotFoundError(
                f'clinicalEnvironmentHoldout.csv not found in productEpi/ or summedEpi/ '
                f'siblings of {pheno_data}. Pass --clinical_source_path explicitly.'
            )

    clinical_cols_present = pd.read_csv(clinical_source_path, nrows=1)\
                              .set_index('IID').columns.tolist()
    df_clin = pd.read_csv(
        os.path.join(results_path, 'participant_environment.csv'),
        usecols=['Participant ID'] + clinical_cols_present
    )
    df_clin.rename(columns={'Participant ID': 'IID'}, inplace=True)

    # Merge
    df = df_prs[['IID', 'PHENOTYPE', prs_col]].merge(df_env, on='IID', how='left')
    df = df.merge(df_clin, on='IID', how='left')

    # Drop rows with any missing in core columns
    core_cols = [prs_col] + [c for c in covariate_cols if c in df.columns] \
                + clinical_cols_present
    n_before = len(df)
    df.dropna(subset=core_cols, inplace=True)
    print(f'Dropped {n_before - len(df)} rows with missing values in core columns.')
    print(f'Final holdout N={len(df)}, cases={df["PHENOTYPE"].sum()}')

    actual_cov = [c for c in covariate_cols if c in df.columns]
    return df, actual_cov, clinical_cols_present


# ---------------------------------------------------------------------------
# results export
# ---------------------------------------------------------------------------

def export_results(results, figPath, pheno):
    rows = []
    for k, v in results.items():
        row = {'model': k}
        row.update(v)
        rows.append(row)
    df_out = pd.DataFrame(rows)
    out = os.path.join(figPath, f'logistic_regression_results.{pheno}.csv')
    df_out.to_csv(out, index=False)
    print(f'Results table saved to {out}')
    return df_out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run_logistic_analysis(pheno, pheno_data, results_path, prs_col,
                          covariate_cols, clinical_source_path=None, n_splits=5):
    figPath = os.path.join(pheno_data, 'figures', 'clinicalFigures')
    os.makedirs(figPath, exist_ok=True)

    df, actual_cov, clinical_cols = load_data(
        pheno_data, results_path, pheno, prs_col,
        covariate_cols, clinical_source_path
    )

    results, base_prob, clin_probs = run_logistic_cv(
        df, prs_col, actual_cov, clinical_cols,
        outcome_col='PHENOTYPE', n_splits=n_splits
    )

    plot_auc_comparison(results, results['_base']['auc'], figPath,
                        prs_col_label=prs_col)
    plot_nri_forest(results, figPath)
    export_results(results, figPath, pheno)


if __name__ == '__main__':
    DEFAULT_COVARIATES = (
        'age,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'
    )

    parser = argparse.ArgumentParser(
        description='Logistic regression: PRSComp + covariates vs + clinical measure'
    )
    parser.add_argument('--pheno_data', help='Path to combinedAnalysis directory')
    parser.add_argument('--pheno',      help='Phenotype name (e.g. type2Diabetes)')
    parser.add_argument('--results_path', help='Path to results directory')
    parser.add_argument('--prs_col',    default='scaled_prs_combined',
                        help='Column name for the composite PRS (default: scaled_prs_combined)')
    parser.add_argument('--covariates', default=DEFAULT_COVARIATES,
                        help='Comma-separated covariate column names')
    parser.add_argument('--clinical_source_path', default=None,
                        help='Path to clinicalEnvironmentHoldout.csv (optional)')
    parser.add_argument('--n_splits',   type=int, default=5,
                        help='Number of CV folds (default: 5)')
    args = parser.parse_args()

    pheno        = args.pheno        or os.environ.get('PHENO')
    pheno_data   = args.pheno_data   or os.environ.get('PHENO_DATA')
    results_path = args.results_path or os.environ.get('RESULTS_PATH')
    prs_col      = args.prs_col      or os.environ.get('PRS_COL', 'scaled_prs_combined')

    if not all([pheno, pheno_data, results_path]):
        raise ValueError(
            'Must provide --pheno, --pheno_data, and --results_path '
            '(or PHENO/PHENO_DATA/RESULTS_PATH env vars).'
        )

    covariate_cols = [c.strip() for c in args.covariates.split(',')]

    run_logistic_analysis(
        pheno, pheno_data, results_path, prs_col,
        covariate_cols, args.clinical_source_path, args.n_splits
    )
