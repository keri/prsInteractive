#!/usr/bin/env python3
# Keri Multerer 2026
# Secondary analysis: combine the five ePRS encodings by fitting them jointly,
# rather than by the max-bin rule used to build ePRS_Comp.
#
# WHY
# ePRS_Comp classifies an individual high risk if ANY retained model places
# them in its top 20%. Five models each flagging about a fifth of the cohort
# gives a union of 64.4% of the holdout -- almost exactly the 67.2% expected if
# the five flags were independent. The rule therefore has a single operating
# point, fixed by that arithmetic rather than chosen, and it cannot express how
# much each encoding contributes.
#
# Fitting one logistic model over the five scores instead returns a calibrated
# probability that can be thresholded anywhere, and gives a coefficient per
# encoding. This is the standard way multi-score PRS are combined and is the
# obvious question a reader will ask of the max-bin rule, so it is reported
# here as a secondary analysis.
#
# VALIDATION DISCIPLINE
# The combining weights are fitted on the validation set and applied once to
# the holdout, which is never used to choose them. Fitting the weights inside
# the holdout by cross-validation gives essentially the same answer (AUC 0.6511
# against 0.6515), so the result is not driven by that choice, but the
# out-of-sample version is the one reported.
#
# AUC CONVENTION
# All AUCs here are computed directly from the score. They are not necessarily
# on the same footing as the AUCs in the main tables, which are produced by
# prsAUCDelongStats.R; the rule-based composite scores 0.605 here against a
# published 0.65, so the two conventions differ and the numbers in this table
# should be compared only with each other.

import argparse
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helper.mcnemar_comparison_calculations import delong_test  # noqa: E402

# Candidate encodings, in a fixed order. Only those retained as statistically
# distinct for a given phenotype are present in the score files -- coeliac
# retains three of the five -- so the set is intersected with what is actually
# there rather than assumed.
CANDIDATE_PRS_COLS = [
    'scaled_prs_main',
    'scaled_prs_epi_summed',
    'scaled_prs_epi_product',
    'scaled_prs_cardio_summed',
    'scaled_prs_cardio_product',
]
PRS_LABELS = {
    'scaled_prs_main':           'PRS_G',
    'scaled_prs_epi_summed':     'PRS_G+G',
    'scaled_prs_epi_product':    'PRS_GxG',
    'scaled_prs_cardio_summed':  'PRS_G+GxE',
    'scaled_prs_cardio_product': 'PRS_GxGxE',
}
FLAG_RATES = [0.644, 0.20, 0.10]


def operating_characteristics(flag, y, base):
    tp = int(((flag == 1) & (y == 1)).sum())
    fp = int(((flag == 1) & (y == 0)).sum())
    fn = int(((flag == 0) & (y == 1)).sum())
    tn = int(((flag == 0) & (y == 0)).sum())
    ppv = tp / max(tp + fp, 1)
    return {
        'flagged': int(flag.sum()), 'flag_rate': float(flag.mean()),
        'sensitivity': tp / max(tp + fn, 1), 'specificity': tn / max(tn + fp, 1),
        'ppv': ppv, 'npv': tn / max(tn + fn, 1),
        'enrichment': ppv / base if base else np.nan,
        'cases_caught': tp,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pheno_data', required=True)
    ap.add_argument('--out_dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    scores = os.path.join(args.pheno_data, 'scores')
    val = pd.read_csv(os.path.join(scores, 'combinedPRSGroups.filtered.csv'))
    hold = pd.read_csv(os.path.join(scores, 'combinedPRSGroups.holdout.filtered.csv'))

    PRS_COLS = [c for c in CANDIDATE_PRS_COLS
                if c in val.columns and c in hold.columns]
    if len(PRS_COLS) < 2:
        raise SystemExit(f'need at least two PRS columns, found {PRS_COLS}')
    dropped = [c for c in CANDIDATE_PRS_COLS if c not in PRS_COLS]
    print(f'encodings retained for this phenotype: '
          f'{[PRS_LABELS[c] for c in PRS_COLS]}')
    if dropped:
        print(f'  not present (not retained as distinct): '
              f'{[PRS_LABELS[c] for c in dropped]}')

    y_val = (val.PHENOTYPE == 2).astype(int).to_numpy()
    y = (hold.PHENOTYPE == 2).astype(int).to_numpy()
    base = y.mean()
    print(f'validation n={len(val)} ({y_val.mean() * 100:.2f}% cases)   '
          f'holdout n={len(hold)} ({base * 100:.2f}% cases)\n')

    # Fit on validation only. Scores are z-scored using validation moments so
    # the holdout contributes nothing to the transform either.
    sc = StandardScaler().fit(val[PRS_COLS])
    X_val = sm.add_constant(sc.transform(val[PRS_COLS]))
    fit = sm.Logit(y_val, X_val).fit(disp=0)
    ens = fit.predict(sm.add_constant(sc.transform(hold[PRS_COLS])))

    # --- coefficients ------------------------------------------------------
    ci = fit.conf_int()
    coef = pd.DataFrame({
        'encoding': ['(intercept)'] + [PRS_LABELS[c] for c in PRS_COLS],
        'coefficient': fit.params,
        'se': fit.bse,
        'ci_low': ci[:, 0], 'ci_high': ci[:, 1],
        'p_value': fit.pvalues,
    })
    coef_path = os.path.join(args.out_dir, 'fitted_prs_ensemble_coefficients.csv')
    coef.to_csv(coef_path, index=False)
    print('Combining weights (validation set, scores z-scored)')
    print(coef.to_string(index=False, float_format=lambda v: f'{v:.4f}'))

    # --- discrimination ----------------------------------------------------
    comparators = {
        'fitted ensemble': ens,
        'ePRS_Comp (max-bin rule)': hold['scaled_prs_combined'].to_numpy(),
        'PRS_G': hold['scaled_prs_main'].to_numpy(),
    }
    print('\nDiscrimination on the untouched holdout')
    for name, s in comparators.items():
        print(f'  {name:<28} AUC {roc_auc_score(y, s):.4f}')

    print()
    for ref in ('ePRS_Comp (max-bin rule)', 'PRS_G'):
        try:
            d = delong_test(y, comparators[ref], ens)
            print(f'  DeLong ensemble vs {ref}: p={d.get("pvalue", np.nan):.3g}')
        except Exception as exc:                                # pragma: no cover
            print(f'  DeLong failed vs {ref}: {exc}')

    # --- operating characteristics at matched flag rates -------------------
    rows = []
    for name, s in comparators.items():
        for rate in FLAG_RATES:
            flag = (s >= np.quantile(s, 1 - rate)).astype(int)
            rows.append({'score': name, 'target_flag_rate': rate,
                         'auc': roc_auc_score(y, s),
                         **operating_characteristics(flag, y, base)})
    # The deployed rule at its own, unchosen operating point.
    rows.append({'score': 'ePRS_Comp binary flag (as deployed)',
                 'target_flag_rate': np.nan,
                 'auc': roc_auc_score(y, hold['combined_high_risk']),
                 **operating_characteristics(
                     hold['combined_high_risk'].to_numpy(), y, base)})
    # Reference: no model at all.
    rows.append({'score': 'flag everyone (no model)', 'target_flag_rate': 1.0,
                 'auc': 0.5,
                 **operating_characteristics(np.ones(len(y), int), y, base)})

    res = pd.DataFrame(rows)
    out = os.path.join(args.out_dir, 'fitted_prs_ensemble_operating.csv')
    res.to_csv(out, index=False)

    print('\n' + '=' * 104)
    print('OPERATING CHARACTERISTICS AT MATCHED FLAG RATES (holdout)')
    print('=' * 104)
    show = ['score', 'flag_rate', 'auc', 'sensitivity', 'specificity',
            'ppv', 'enrichment', 'cases_caught']
    print(res[show].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'\n-> {coef_path}\n-> {out}')


if __name__ == '__main__':
    main()
