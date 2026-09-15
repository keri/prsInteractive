#!/usr/bin/env python3
# Keri Multerer 2026
# Nested reclassification of a clinical risk model with and without ePRS_Comp.
#
# WHY THIS EXISTS
# The published Figure 6 NRI compares a clinical threshold rule (e.g. glucose
# below a cut-point) against a PRS rule (combined_high_risk). Those are two
# unrelated classifiers, not nested models, and they flag very different
# fractions of the cohort -- the glucose rule flags under 1%, combined_high_risk
# flags 64.4% of the holdout. With a budget mismatch that large the event and
# non-event components are arithmetically forced: almost every case moves up
# and so does almost every control. Pencina's NRI is defined for nested models,
# so the comparison has to be reference model versus reference model plus the
# new marker.
#
# WHAT IT DOES
# Fits a clinical model and the same model plus ePRS_Comp, pools out-of-fold
# predicted probabilities over 5-fold stratified CV, and reclassifies between
# them. Two structures are produced because they answer different questions:
#
#   combined     age + sex + HbA1c + BMI + glucose, versus that plus ePRS.
#                The realistic clinical comparator: what does ePRS add to
#                everything a clinician already has?
#   per-measure  age + sex + one clinical measure, versus that plus ePRS,
#                repeated for each measure. Shows which single measure ePRS
#                complements, at the cost of three tests instead of one.
#
# ANCESTRY PCs are excluded from both models. They are not clinical, and
# putting them in only the ePRS arm would let the comparison credit ePRS with
# ancestry adjustment. Both arms therefore carry age and sex alone as
# non-clinical terms, so the only difference between them is ePRS.
#
# REPORTED QUANTITIES
#   NRI at a fixed 20% classification budget, with the event and non-event
#       components separated. Both arms are binarised at the SAME threshold --
#       the reference model's 80th percentile -- so the budget is set once and
#       ePRS has to earn reclassification rather than being handed a larger
#       one. 20% matches the risk stratum used throughout the paper. The value
#       at the Youden point is reported alongside for reference.
#   Category-free NRI, which needs no threshold at all: the net proportion of
#       cases whose predicted risk rises, plus the net proportion of controls
#       whose predicted risk falls. Reported because any threshold choice is
#       arguable and cfNRI removes that degree of freedom.
#   IDI, the mean change in predicted risk separation, which is the standard
#       companion to NRI and is not threshold-dependent.
#   AUC for both arms with a DeLong test on the difference.
#
# GLUCOSE THRESHOLD
# The binary clinical flags used by the old analysis took glucose at
# 12 mmol/L. That is above the 99th percentile of the holdout and above the
# diagnostic cut-point on any reading, so it classified 99.4% of measured
# participants -- including 2,206 of 2,696 cases -- as clinically normal.
# UK Biobank glucose is non-fasting, so the diagnostic threshold is
# 11.1 mmol/L. This script models the continuous measures directly and so does
# not depend on a cut-point, but it reports the flags under both values for
# reference.

import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helper.mcnemar_comparison_calculations import delong_test  # noqa: E402

CLINICAL = {
    'Glycated haemoglobin (HbA1c)': 'hba1c',
    'Body mass index (BMI)':        'bmi',
    'Glucose':                      'glucose',
}

# Published diagnostic thresholds, for the reference table only. Glucose is
# the UK Biobank non-fasting measurement, so 11.1 mmol/L is the diagnostic
# cut-point; 12 is the value the previous analysis used.
THRESHOLDS = {
    'Glycated haemoglobin (HbA1c)': [('prediabetes', 42)],
    'Body mass index (BMI)':        [('overweight', 25)],
    'Glucose':                      [('diabetic (random)', 11.1), ('previously used', 12)],
}

COVARIATES = ['age', 'SEX']


# ---------------------------------------------------------------------------
# reclassification statistics
# ---------------------------------------------------------------------------

def nri_categorical(y, prob_ref, prob_new, threshold):
    """Pencina NRI with both arms binarised at the reference model's threshold."""
    case = y.astype(bool)
    old = prob_ref >= threshold
    new = prob_new >= threshold

    up_c   = int((new & ~old & case).sum())
    down_c = int((~new & old & case).sum())
    up_n   = int((new & ~old & ~case).sum())
    down_n = int((~new & old & ~case).sum())
    n_c, n_n = int(case.sum()), int((~case).sum())

    ev = (up_c - down_c) / n_c
    ne = (down_n - up_n) / n_n
    se_ev = np.sqrt(up_c + down_c) / n_c
    se_ne = np.sqrt(up_n + down_n) / n_n
    se = np.sqrt(se_ev ** 2 + se_ne ** 2)
    return {
        'nri': ev + ne, 'nri_events': ev, 'nri_non_events': ne,
        'ci_low': ev + ne - 1.96 * se, 'ci_high': ev + ne + 1.96 * se,
        'se': se, 'threshold': threshold,
        'up_cases': up_c, 'down_cases': down_c,
        'up_controls': up_n, 'down_controls': down_n,
        'n_flagged_ref': int(old.sum()), 'n_flagged_new': int(new.sum()),
    }


def nri_category_free(y, prob_ref, prob_new):
    """Category-free NRI: needs no threshold, so no cut-point to argue over."""
    case = y.astype(bool)
    d = prob_new - prob_ref
    ev = (np.sum(d[case] > 0) - np.sum(d[case] < 0)) / case.sum()
    ne = (np.sum(d[~case] < 0) - np.sum(d[~case] > 0)) / (~case).sum()
    se = np.sqrt(1 / case.sum() + 1 / (~case).sum())
    return {'cfnri': ev + ne, 'cfnri_events': ev, 'cfnri_non_events': ne,
            'cfnri_ci_low': ev + ne - 1.96 * se,
            'cfnri_ci_high': ev + ne + 1.96 * se}


def integrated_discrimination(y, prob_ref, prob_new):
    """IDI: gain in mean predicted risk separation between cases and controls."""
    case = y.astype(bool)
    idi = ((prob_new[case].mean() - prob_ref[case].mean())
           - (prob_new[~case].mean() - prob_ref[~case].mean()))
    return {'idi': idi}


def youden_threshold(y, prob):
    fpr, tpr, thr = roc_curve(y, prob)
    return float(thr[np.argmax(tpr - fpr)])


# ---------------------------------------------------------------------------
# model fitting
# ---------------------------------------------------------------------------

def pooled_probabilities(df, features, y, seed=42, n_splits=5):
    """Out-of-fold predicted probabilities, features z-scored within each fold."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    cont = [c for c in features if c != 'SEX']
    out = np.full(len(df), np.nan)

    for tr, te in skf.split(np.arange(len(df)), y):
        d_tr, d_te = df.iloc[tr].copy(), df.iloc[te].copy()
        if cont:
            sc = StandardScaler()
            d_tr[cont] = sc.fit_transform(d_tr[cont])
            d_te[cont] = sc.transform(d_te[cont])
        m = LogisticRegression(penalty=None, solver='lbfgs', max_iter=2000,
                               random_state=seed)
        m.fit(d_tr[features].values, y[tr])
        out[te] = m.predict_proba(d_te[features].values)[:, 1]
    return out


def compare(df, y, clinical_features, prs_col, label, seed, flag_rate=0.20):
    """Clinical model versus the same model plus ePRS.

    Both arms are binarised at the SAME threshold, taken from the reference
    model, so the classification budget is set once and ePRS has to earn
    reclassification rather than being handed a larger one.

    The threshold is the reference model's `flag_rate` quantile -- the top 20%
    by default, matching the risk stratum used throughout the paper -- rather
    than its Youden point. Youden weights a missed case and a false alarm
    equally, which is a different and unstated assumption, and it flags a
    different fraction of the cohort for every clinical model, so the budget
    would vary between rows of the results table. It is still reported
    alongside for reference.
    """
    ref_feats = COVARIATES + clinical_features
    new_feats = ref_feats + [prs_col]

    p_ref = pooled_probabilities(df, ref_feats, y, seed)
    p_new = pooled_probabilities(df, new_feats, y, seed)

    thr = float(np.quantile(p_ref, 1 - flag_rate))
    thr_youden = youden_threshold(y, p_ref)
    auc_ref = roc_auc_score(y, p_ref)
    auc_new = roc_auc_score(y, p_new)

    row = {'model': label,
           'clinical_terms': ' + '.join(clinical_features),
           'flag_rate': flag_rate,
           'auc_clinical': auc_ref, 'auc_clinical_plus_eprs': auc_new,
           'delta_auc': auc_new - auc_ref}
    row.update(nri_categorical(y, p_ref, p_new, thr))
    yd = nri_categorical(y, p_ref, p_new, thr_youden)
    row['nri_at_youden'] = yd['nri']
    row['youden_flag_rate'] = yd['n_flagged_ref'] / len(y)
    row.update(nri_category_free(y, p_ref, p_new))
    row.update(integrated_discrimination(y, p_ref, p_new))

    # DeLong on the two correlated AUCs; the helper prints its own report.
    try:
        d = delong_test(y, p_ref, p_new)
        row['delong_p'] = d.get('pvalue', np.nan)
    except Exception as exc:                                   # pragma: no cover
        print(f'  [warn] DeLong failed for {label}: {exc}')
        row['delong_p'] = np.nan
    return row


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pheno_data', required=True,
                    help='.../results/<pheno>/combinedAnalysis')
    ap.add_argument('--results_path', required=True,
                    help='directory holding participant_environment.csv')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--prs_col', default='scaled_prs_ensemble',
                    help='ePRS_Comp is the fitted ensemble, so its score lives '
                         'in --extra_scores_file rather than in the holdout '
                         'file. Pass bin_combined to use the rule instead.')
    ap.add_argument('--extra_scores_file', default=None,
                    help='fitted_ensemble_holdout_scores.csv, merged on IID.')
    ap.add_argument('--covar_file', default=None,
                    help='Overrides the default holdout covariate PRS file.')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--flag_rate', type=float, default=0.20,
                    help='Fraction of the cohort the reference model flags; '
                         'both arms are binarised at that same threshold. '
                         'Default 0.20 to match the risk stratum used '
                         'throughout the paper.')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    holdout = pd.read_csv(
        os.path.join(args.pheno_data, 'scores',
                     'combinedPRSGroups.holdout.filtered.csv'))
    if args.extra_scores_file:
        extra = pd.read_csv(args.extra_scores_file)
        cols = ['IID'] + [c for c in extra.columns
                          if c != 'IID' and c not in holdout.columns]
        holdout = holdout.merge(extra[cols], on='IID', how='left')
    env = pd.read_csv(os.path.join(args.results_path, 'participant_environment.csv'),
                      usecols=['Participant ID'] + list(CLINICAL))
    env = env.rename(columns={'Participant ID': 'IID'})

    # age/SEX come from the holdout covariate PRS file, which is the same
    # source the rest of the pipeline uses for these terms.
    covar_path = args.covar_file or os.path.join(
        args.pheno_data, 'scores', 'prsScores', 'covariate.holdout.mixed.prs.csv')
    keep = ['IID', 'PHENOTYPE', args.prs_col]
    df = holdout[[c for c in keep if c in holdout.columns]].merge(env, on='IID', how='left')
    if not os.path.exists(covar_path):
        raise SystemExit(f'covariate file not found: {covar_path}')
    cov = pd.read_csv(covar_path)
    take = ['IID'] + [c for c in COVARIATES if c in cov.columns]
    df = df.merge(cov[take].drop_duplicates('IID'), on='IID', how='left')
    missing_cov = [c for c in COVARIATES if c not in df.columns]
    if missing_cov:
        raise SystemExit(f'covariates not found: {missing_cov} (looked in {covar_path})')

    df['PHENOTYPE'] = df['PHENOTYPE'].map({1: 0, 2: 1}).astype(float)

    # Complete cases only, so both arms of every comparison are fitted on the
    # same individuals and the AUCs are directly comparable.
    need = ['PHENOTYPE', args.prs_col] + COVARIATES + list(CLINICAL)
    before = len(df)
    df = df.dropna(subset=need).reset_index(drop=True)
    y = df['PHENOTYPE'].to_numpy()
    print(f'holdout {before} -> {len(df)} complete cases '
          f'({int(y.sum())} cases, {int((1 - y).sum())} controls)\n')

    print('Clinical threshold reference (continuous measures are modelled '
          'directly; these are for context only)')
    for col, flags in THRESHOLDS.items():
        for name, t in flags:
            below = df[col] <= t
            cases_below = int(((y == 1) & below).sum())
            print(f'  {col} <= {t:<5} [{name}]: {below.sum():6d} '
                  f'({below.mean() * 100:5.1f}%) | cases below {cases_below}')
    print()

    rows = []
    rows.append(compare(df, y, list(CLINICAL), args.prs_col,
                        'combined clinical model', args.seed, args.flag_rate))
    for col in CLINICAL:
        rows.append(compare(df, y, [col], args.prs_col,
                            f'per-measure: {CLINICAL[col]}', args.seed,
                            args.flag_rate))

    res = pd.DataFrame(rows)
    out = os.path.join(args.out_dir, 'nested_nri_clinical_eprs.csv')
    res.to_csv(out, index=False)

    print('\n' + '=' * 78)
    print('NESTED RECLASSIFICATION: clinical model vs clinical model + ePRS_Comp')
    print('=' * 78)
    show = ['model', 'auc_clinical', 'auc_clinical_plus_eprs', 'delta_auc',
            'delong_p', 'nri', 'nri_events', 'nri_non_events', 'cfnri', 'idi',
            'nri_at_youden']
    print(res[show].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'\n-> {out}')


if __name__ == '__main__':
    main()
