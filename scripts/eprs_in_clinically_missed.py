#!/usr/bin/env python3
# Keri Multerer 2026
# What does ePRS_Comp contribute among people clinical measures do not flag?
#
# WHY THIS IS A SEPARATE QUESTION
# The nested reclassification analysis asks whether ePRS improves a clinical
# model across everyone. It came out near null, but that is dominated by the
# people clinical measures already identify: HbA1c and glucose are diagnostic
# criteria for T2D measured at the same baseline visit as the prevalent-disease
# outcome, so among biomarker-positive individuals the clinical model is close
# to re-deriving the case definition and nothing can add to it.
#
# The clinically useful question is conditional: among people the clinical
# measures do NOT flag, does ePRS still separate cases from controls? That is
# not answerable from an overall NRI, because the overall figure is an average
# dominated by the subgroup where clinical information is near-perfect.
#
# WHY NOT THE OLD FIGURE 6 RIGHT PANEL
# It attempted this subgroup analysis, but defined "clinically low risk" with
# glucose <= 12 mmol/L, which is above the diagnostic cut-point and left 99.4%
# of participants in the "low risk" group -- so the subgroup was the whole
# cohort, not the missed population. It then compared a threshold rule flagging
# under 1% against a PRS rule flagging 64.4%, so the reclassification counts
# were set by the budget mismatch rather than by discrimination.
#
# HOW THE SUBGROUP IS DEFINED
# Five definitions are reported. Each measure is shown on its own because they
# have very different blind spots, and two joint definitions bracket what
# "clinically missed" can mean:
#
#   HbA1c-negative       HbA1c < 48 mmol/mol, the diagnostic cut-point.
#   glucose-negative     glucose < 11.1 mmol/L. Barely a restriction here: the
#                        UK Biobank measurement is non-fasting and rarely
#                        reaches the diagnostic value, so this keeps 99% of
#                        participants.
#   BMI-negative         BMI < 25. Not a diagnostic criterion but a risk
#                        factor, so this is the group where clinical suspicion
#                        is lowest rather than where testing was negative.
#   biomarker-negative   both diagnostic cut-points together: people whose
#                        routine diabetes bloods would not have diagnosed them.
#   model-negative       below the operating threshold of a fitted clinical
#                        model (age + sex + HbA1c + BMI + glucose). Uses the
#                        measures jointly and continuously, so someone with
#                        sub-diagnostic bloods plus high BMI and older age is
#                        excluded even though every cut-point calls them
#                        negative. This is the stricter and more realistic
#                        test, and it is close to a subset of
#                        biomarker-negative; the overlap table prints how the
#                        two differ and how many cases sit in the gap.
#
# WHAT IS MEASURED, AND AT A MATCHED BUDGET
# Discrimination within the subgroup (AUC, with a DeLong test of ePRS_Comp
# against PRS_G, since the paper's claim is that the composite beats the
# additive score). Then, at a FIXED flag rate inside the subgroup -- top 10%
# and top 20% -- the number of cases captured, the positive predictive value,
# and the enrichment over the subgroup's own base rate. Fixing the flag rate is
# the point: it is what stops a score from looking good merely by flagging more
# people, which is what went wrong in the published analysis. Note that the
# comparison BETWEEN scores is insensitive to that rate, since both are held to
# it; the subgroup boundary, by contrast, is a real assumption that changes who
# is analysed, which is why several are reported rather than one.
#
# Ranking uses the CONTINUOUS composite score, not combined_high_risk. The
# binary flag marks 64.4% of the holdout, so it cannot express a top decile.

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


# ---------------------------------------------------------------------------
# clinical measure discovery
# ---------------------------------------------------------------------------

CLIN_DIR = 'clinicalMeasures'
CLIN_CSV = 'combinedPRS.holdout.ClinicalMeasures.csv'
CLIN_XLSX = 'clinicalMeasurePerformanceResultsholdout.bin.xlsx'


def discover_clinical(pheno_data):
    """Find this phenotype's clinical measures, thresholds and risk directions.

    The measures differ by disease -- HbA1c, BMI and glucose for T2D; basal
    metabolic rate, urea and haemoglobin for coeliac -- and the direction that
    counts as high risk differs with them: T2D measures are elevated in cases,
    the coeliac measures are reduced, consistent with malabsorption. Hard-coding
    either set makes the analysis silently wrong for the other disease, so both
    are read from the pipeline's own outputs.

    Returns (measures, thresholds, directions) where directions[m] is True when
    HIGH values indicate risk.
    """
    base = os.path.join(pheno_data, 'scores', CLIN_DIR)
    csv_path = os.path.join(base, CLIN_CSV)
    if not os.path.exists(csv_path):
        raise SystemExit(f'clinical measures file not found: {csv_path}')

    cols = pd.read_csv(csv_path, nrows=0).columns
    measures = [c[:-len('_binary')] for c in cols
                if c.endswith('_binary') and c[:-len('_binary')] in cols]
    if not measures:
        raise SystemExit(f'no paired <measure>/<measure>_binary columns in {csv_path}')

    thresholds, directions = {}, {}
    xlsx = os.path.join(base, CLIN_XLSX)
    if os.path.exists(xlsx):
        t = pd.read_excel(xlsx, sheet_name='thresholds RiskDirection')
        for _, r in t.iterrows():
            m = r['clinical_measure']
            thresholds[m] = float(r['thresholds_used'])
            directions[m] = bool(r['risk_directions'])
    missing = [m for m in measures if m not in thresholds]
    if missing:
        raise SystemExit(f'no threshold recorded for {missing} in {xlsx}')

    print(f'clinical measures for this phenotype ({len(measures)}):')
    for m in measures:
        if directions[m]:
            print(f'  {m:<30} high values indicate risk, '
                  f'threshold {thresholds[m]:g}')
        else:
            print(f'  {m:<30} low values indicate risk, so the high-risk '
                  f'group is the bottom {LOW_RISK_QUANTILE:.0%} '
                  f'(stored threshold {thresholds[m]:g} sits on the wrong tail)')
    return measures, thresholds, directions


LOW_RISK_QUANTILE = 0.20


def clinical_high_risk(df, measure, thresholds, directions,
                       q=LOW_RISK_QUANTILE):
    """Boolean flag: is this individual high risk on this measure?

    Where HIGH values indicate risk the stored threshold is used directly:
    for T2D these are the published clinical cut-points (HbA1c 41 mmol/mol,
    BMI 25, glucose 12) and need no adjustment.

    Where LOW values indicate risk the stored threshold cannot be used. Those
    thresholds are derived rather than supplied, and convert_to_binary places
    them at the 80th percentile irrespective of direction -- so for coeliac,
    where reduced basal metabolic rate, urea and haemoglobin indicate
    malabsorption, "below the threshold" captures the bottom 80% of the cohort
    rather than the bottom 20%. That flags four-fifths of everyone as high
    risk at a case rate of 1.18-1.22% against a 1.11% baseline, which is
    almost no enrichment, and leaves the complementary "clinically missed"
    group too small to analyse. Taking the bottom 20% instead gives case rates
    of 1.63-1.84%, a 1.5- to 1.7-fold enrichment.

    The stored value is therefore replaced by the qth percentile when the
    direction is low, and the substitution is reported.
    """
    if directions[measure]:
        return df[measure] > thresholds[measure]
    return df[measure] < df[measure].quantile(q)


COVARIATES = ['age', 'SEX']
FLAG_RATES = [0.10, 0.20]


def odds_ratio_ci(flagged, y):
    """OR for case status, flagged vs not, with a Woolf 95% CI."""
    a = int(((flagged == 1) & (y == 1)).sum())
    b = int(((flagged == 1) & (y == 0)).sum())
    c = int(((flagged == 0) & (y == 1)).sum())
    d = int(((flagged == 0) & (y == 0)).sum())
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    or_ = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return or_, or_ * np.exp(-1.96 * se), or_ * np.exp(1.96 * se)


def clinical_model_risk(df, y, measures, seed=42, n_splits=5):
    """Out-of-fold risk from age + sex + this phenotype's clinical measures."""
    feats = COVARIATES + list(measures)
    cont = [c for c in feats if c != 'SEX']
    out = np.full(len(df), np.nan)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in skf.split(np.arange(len(df)), y):
        a, b = df.iloc[tr].copy(), df.iloc[te].copy()
        sc = StandardScaler()
        a[cont] = sc.fit_transform(a[cont])
        b[cont] = sc.transform(b[cont])
        m = LogisticRegression(penalty=None, solver='lbfgs', max_iter=2000,
                               random_state=seed)
        m.fit(a[feats].values, y[tr])
        out[te] = m.predict_proba(b[feats].values)[:, 1]
    return out


def evaluate(sub, y, scores, label, rows):
    """Discrimination and fixed-budget yield for each score in the subgroup."""
    base = y.mean()
    for name, s in scores.items():
        auc = roc_auc_score(y, s)
        row = {'subgroup': label, 'score': name, 'n': len(y),
               'n_cases': int(y.sum()), 'base_rate': base, 'auc': auc}
        for rate in FLAG_RATES:
            cut = np.quantile(s, 1 - rate)
            flag = (s >= cut).astype(int)
            caught = int(((flag == 1) & (y == 1)).sum())
            ppv = caught / max(flag.sum(), 1)
            or_, lo, hi = odds_ratio_ci(flag, y)
            tag = f'top{int(rate * 100)}'
            row[f'{tag}_n_flagged'] = int(flag.sum())
            row[f'{tag}_cases_caught'] = caught
            row[f'{tag}_ppv'] = ppv
            row[f'{tag}_enrichment'] = ppv / base if base else np.nan
            row[f'{tag}_or'] = or_
            row[f'{tag}_or_lo'] = lo
            row[f'{tag}_or_hi'] = hi
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pheno_data', required=True)
    ap.add_argument('--results_path', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--composite_col', default='scaled_prs_ensemble',
                    help='ePRS_Comp is the fitted ensemble; pass '
                         'scaled_prs_combined to use the rule instead.')
    ap.add_argument('--extra_scores_file', default=None,
                    help='fitted_ensemble_holdout_scores.csv, merged on IID.')
    ap.add_argument('--additive_col', default='scaled_prs_main')
    ap.add_argument('--covar_file', default=None)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    h = pd.read_csv(os.path.join(args.pheno_data, 'scores',
                                 'combinedPRSGroups.holdout.filtered.csv'))
    if args.extra_scores_file:
        extra = pd.read_csv(args.extra_scores_file)
        cols = ['IID'] + [c for c in extra.columns
                          if c != 'IID' and c not in h.columns]
        h = h.merge(extra[cols], on='IID', how='left')
    measures, thresholds, directions = discover_clinical(args.pheno_data)
    env = pd.read_csv(os.path.join(args.pheno_data, 'scores', CLIN_DIR, CLIN_CSV),
                      usecols=['IID'] + measures)
    covar_path = args.covar_file or os.path.join(
        args.pheno_data, 'scores', 'prsScores', 'covariate.holdout.mixed.prs.csv')
    cov = pd.read_csv(covar_path)

    keep = ['IID', 'PHENOTYPE', args.composite_col, args.additive_col]
    df = (h[[c for c in keep if c in h.columns]]
          .merge(env, on='IID', how='left')
          .merge(cov[['IID'] + [c for c in COVARIATES if c in cov.columns]]
                 .drop_duplicates('IID'), on='IID', how='left'))

    need = ['PHENOTYPE', args.composite_col, args.additive_col] + measures + COVARIATES
    before = len(df)
    df = df.dropna(subset=need).reset_index(drop=True)
    y = (df.PHENOTYPE == 2).astype(int).to_numpy()
    print(f'holdout {before} -> {len(df)} complete cases '
          f'({int(y.sum())} cases, base rate {y.mean() * 100:.2f}%)\n')

    # --- subgroup definitions -------------------------------------------------
    # Not flagged by each measure individually, and by all of them together.
    per_measure = {m: ~clinical_high_risk(df, m, thresholds, directions).to_numpy()
                   for m in measures}
    biomarker_neg = np.logical_and.reduce(list(per_measure.values()))

    clin_risk = clinical_model_risk(df, y, measures, args.seed)
    fpr, tpr, thr = roc_curve(y, clin_risk)
    clin_thr = float(thr[np.argmax(tpr - fpr)])
    model_neg = clin_risk < clin_thr

    # Each biomarker is reported on its own as well as jointly, because they
    # have very different blind spots: HbA1c alone is the strongest single
    # discriminator, glucose is non-fasting and so rarely diagnostic here, and
    # BMI is a risk factor rather than a diagnostic criterion.
    subgroups = ([(f'not flagged: {m}', v) for m, v in per_measure.items()]
                 + [('biomarker-negative', biomarker_neg),
                    ('model-negative', model_neg)])

    print('Subgroups (individuals clinical measures would not flag)')
    for lab, mask in subgroups:
        n, nc = int(mask.sum()), int(y[mask].sum())
        print(f'  {lab:<20} n={n:6d} ({mask.mean() * 100:4.1f}%)  '
              f'cases={nc:5d} ({nc / n * 100:5.2f}%)  '
              f'cases missed by clinical = {nc} of {int(y.sum())}')
    print()

    rows = []
    scores_all = {'ePRS_Comp': df[args.composite_col].to_numpy(),
                  'PRS_G': df[args.additive_col].to_numpy()}
    rows = evaluate(None, y, scores_all, 'whole holdout', rows)
    print('\nOverlap of the two headline definitions:')
    bn, mn = biomarker_neg, model_neg
    print(f'  both negative           {int((bn & mn).sum()):6d}  '
          f'cases {int(y[bn & mn].sum()):5d}')
    print(f'  biomarker-neg only      {int((bn & ~mn).sum()):6d}  '
          f'cases {int(y[bn & ~mn].sum()):5d}   '
          f'(sub-diagnostic values plus age/BMI push model risk up)')
    print(f'  model-neg only          {int((~bn & mn).sum()):6d}  '
          f'cases {int(y[~bn & mn].sum()):5d}')
    print(f'  both positive           {int((~bn & ~mn).sum()):6d}  '
          f'cases {int(y[~bn & ~mn].sum()):5d}\n')

    for lab, mask in subgroups:
        sc = {k: v[mask] for k, v in scores_all.items()}
        rows = evaluate(None, y[mask], sc, lab, rows)

        # Does the composite beat the additive score where it matters?
        try:
            d = delong_test(y[mask], sc['PRS_G'], sc['ePRS_Comp'])
            print(f'  [{lab}] DeLong ePRS_Comp vs PRS_G: '
                  f'p={d.get("pvalue", np.nan):.3g}')
        except Exception as exc:                                # pragma: no cover
            print(f'  [{lab}] DeLong failed: {exc}')

    res = pd.DataFrame(rows)
    out = os.path.join(args.out_dir, 'eprs_in_clinically_missed.csv')
    res.to_csv(out, index=False)

    print('\n' + '=' * 96)
    print('ePRS_Comp AMONG INDIVIDUALS CLINICAL MEASURES DO NOT FLAG')
    print('(fixed flag rate inside each subgroup, so scores are compared at '
          'the same classification budget)')
    print('=' * 96)
    show = ['subgroup', 'score', 'n', 'n_cases', 'base_rate', 'auc',
            'top20_cases_caught', 'top20_ppv', 'top20_enrichment',
            'top20_or', 'top20_or_lo', 'top20_or_hi']
    print(res[show].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'\n-> {out}')


if __name__ == '__main__':
    main()
