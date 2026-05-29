#!/usr/bin/env python3
"""
summarize_env_by_prs_group.py
-----------------------------
Summarize environmental/clinical feature distributions across PRS risk groups,
for use in addressing coauthor questions about E-feature timing and circularity.

For each PRS model in combinedPRSGroups.csv (and .holdout.csv separately):

  Groups
  ------
  controls:        PHENOTYPE == 1
  all_cases:       PHENOTYPE == 2
  high_risk_cases: PHENOTYPE == 2  AND  top-decile PRS bin (bin == 10)

  Metrics per group
  -----------------
  n, pct_male, age_at_participation_mean
  age_E11_diagnosis_mean  (year of E11 first report - year of birth)
  age_K90_diagnosis_mean  (year of K90 first report - year of birth)
    NOTE: participant.csv contains K90 (intestinal malabsorption), not K10.
    UK Biobank measures E features at baseline recruitment — age_at_participation
    and age_at_measurement are therefore the same timepoint.
  mean of each E (clinical/environmental) feature used by that model

  E features are derived from featureScoresReducedFinalModel.combined.filtered.csv:
  the leading token of cardio interaction features (e.g. "Diastolic blood pressure"
  from "Diastolic blood pressure,rs1234,rs5678"). Intercept, age, SEX, and PC1-10
  are excluded.

Output
------
  <pheno_path>/env_summary_by_prs_group.csv
  <pheno_path>/env_summary_by_prs_group.holdout.csv
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

SKIP_TERMS = {'(Intercept)', 'age', 'SEX'} | {f'PC{i}' for i in range(1, 11)}

BIN_TO_MODEL = {
    'bin_epi_product':    'epi_product',
    'bin_cardio_product': 'cardio_product',
    'bin_main':           'main',
    'bin_epi_summed':     'epi_summed',
    'bin_cardio_summed':  'cardio_summed',
}

SCALED_TO_BIN = {
    'scaled_prs_epi_product':    'bin_epi_product',
    'scaled_prs_cardio_product': 'bin_cardio_product',
    'scaled_prs_main':           'bin_main',
    'scaled_prs_epi_summed':     'bin_epi_summed',
    'scaled_prs_cardio_summed':  'bin_cardio_summed',
}


def load_demographics(data_path: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(data_path, 'participant.csv'))
    df = df.rename(columns={'Participant ID': 'IID'})

    for code, col in [
        ('E11', 'Date E11 first reported (non-insulin-dependent diabetes mellitus)'),
        ('K90', 'Date K90 first reported (intestinal malabsorption)'),
    ]:
        if col in df.columns:
            diag_year = pd.to_datetime(df[col], errors='coerce').dt.year
            df[f'age_{code}_diagnosis'] = diag_year - df['Year of birth']

    df['is_male'] = (df['Sex'].str.strip().str.lower() == 'male').astype(float)
    df = df.rename(columns={'Age at recruitment': 'age_at_participation'})

    keep = ['IID', 'is_male', 'age_at_participation',
            'age_E11_diagnosis', 'age_K90_diagnosis']
    return df[[c for c in keep if c in df.columns]]


_RS_RE = re.compile(r'^rs\d+$', re.IGNORECASE)


def get_model_e_features(combined_path: str) -> dict[str, list[str]]:
    """Return {model: [E feature names]} from the combined model file.

    Only returns the leading clinical-variable token of cardio interaction
    features (e.g. 'Diastolic blood pressure' from
    'Diastolic blood pressure,rs1234,rs5678').  rsID-only pair features
    (epi) and covariate terms are excluded.
    """
    combined = pd.read_csv(combined_path)
    model_e: dict[str, list[str]] = {}
    for model in combined['model'].unique():
        sub = combined[combined['model'] == model]
        e_vars = sorted({
            tok
            for f in sub['feature'].astype(str)
            if ',' in f
            for tok in [f.split(',')[0]]
            if tok not in SKIP_TERMS and not _RS_RE.match(tok)
        })
        model_e[model] = e_vars
    return model_e


def add_bins_from_scaled(df: pd.DataFrame) -> pd.DataFrame:
    """Compute decile bins from scaled_prs_* when bin_* columns are absent."""
    df = df.copy()
    for scaled_col, bin_col in SCALED_TO_BIN.items():
        if scaled_col in df.columns and bin_col not in df.columns:
            df[bin_col] = pd.qcut(
                df[scaled_col], q=10,
                labels=range(1, 11), duplicates='drop',
            ).astype(float)
    return df


def _mean(series: pd.Series) -> float | None:
    v = series.dropna()
    return round(float(v.mean()), 2) if len(v) > 0 else None


def _pct_male(series: pd.Series) -> float | None:
    v = series.dropna()
    return round(float(v.mean() * 100), 1) if len(v) > 0 else None


def summarise(df: pd.DataFrame, model_e_features: dict, env_cols: list) -> pd.DataFrame:
    records = []

    for bin_col, model in BIN_TO_MODEL.items():
        if bin_col not in df.columns:
            continue

        e_features = [f for f in model_e_features.get(model, []) if f in env_cols]

        groups = {
            'controls':        df[df['PHENOTYPE'] == 1],
            'all_cases':       df[df['PHENOTYPE'] == 2],
            'high_risk_cases': df[(df['PHENOTYPE'] == 2) & (df[bin_col] == 10)],
        }

        for group_name, grp in groups.items():
            rec: dict = {'model': model, 'group': group_name, 'n': len(grp)}

            rec['pct_male']                  = _pct_male(grp.get('is_male', pd.Series(dtype=float)))
            rec['age_at_participation_mean'] = _mean(grp.get('age_at_participation', pd.Series(dtype=float)))
            rec['age_E11_diagnosis_mean']    = _mean(grp.get('age_E11_diagnosis', pd.Series(dtype=float)))
            rec['age_K90_diagnosis_mean']    = _mean(grp.get('age_K90_diagnosis', pd.Series(dtype=float)))

            for feat in e_features:
                if feat in grp.columns:
                    rec[f'{feat}_mean'] = _mean(grp[feat])

            records.append(rec)

    return pd.DataFrame(records)


def process_file(
    groups_path: str,
    demog: pd.DataFrame,
    env: pd.DataFrame,
    model_e_features: dict,
) -> pd.DataFrame:
    df = pd.read_csv(groups_path)
    df = add_bins_from_scaled(df)

    env_renamed = env.rename(columns={'Participant ID': 'IID'})
    df = df.merge(demog,        on='IID', how='left')
    df = df.merge(env_renamed,  on='IID', how='left')

    env_cols = [c for c in env.columns if c != 'Participant ID']
    return summarise(df, model_e_features, env_cols)


def main(pheno_path: str, data_path: str, combined_path: str) -> None:
    print(f"Loading demographics : {data_path}/participant.csv")
    demog = load_demographics(data_path)
    print(f"  {len(demog)} participants")

    env_path = os.path.join(data_path, 'participant_environment.csv')
    print(f"Loading environment  : {env_path}")
    env = pd.read_csv(env_path)
    print(f"  {len(env)} participants, {len(env.columns) - 1} E features")

    print(f"Loading combined model: {combined_path}")
    model_e_features = get_model_e_features(combined_path)
    for m, feats in sorted(model_e_features.items()):
        if feats:
            print(f"  {m}: {len(feats)} E features")

    scores_dir = os.path.join(pheno_path, 'combinedAnalysis', 'scores')

    for fname, label, out_name in [
        ('combinedPRSGroups.csv',         'discovery', 'env_summary_by_prs_group.csv'),
        ('combinedPRSGroups.holdout.csv', 'holdout',   'env_summary_by_prs_group.holdout.csv'),
    ]:
        fpath = os.path.join(scores_dir, fname)
        if not os.path.exists(fpath):
            print(f"\n[skip] not found: {fpath}")
            continue

        print(f"\nProcessing {label}: {fpath}")
        summary = process_file(fpath, demog, env, model_e_features)

        out_path = os.path.join(pheno_path, out_name)
        summary.to_csv(out_path, index=False)
        print(f"  Wrote {len(summary)} rows → {out_path}")

        preview_cols = ['model', 'group', 'n', 'pct_male',
                        'age_at_participation_mean',
                        'age_E11_diagnosis_mean', 'age_K90_diagnosis_mean']
        preview_cols = [c for c in preview_cols if c in summary.columns]
        print(summary[preview_cols].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarise E features by PRS risk group for coauthor review"
    )
    parser.add_argument(
        '--pheno_path', required=True,
        help='Path to phenotype results directory (e.g. results/type2Diabetes)',
    )
    parser.add_argument(
        '--data_path',
        default='/Users/kerimulterer/prsInteractive/data',
        help='Path to data directory containing participant.csv and participant_environment.csv',
    )
    parser.add_argument(
        '--combined_path',
        default=None,
        help='Path to featureScoresReducedFinalModel.combined.filtered.csv '
             '(default: <pheno_path>/combinedAnalysis/scores/...)',
    )
    args = parser.parse_args()

    combined_path = args.combined_path or os.path.join(
        args.pheno_path, 'combinedAnalysis', 'scores',
        'featureScoresReducedFinalModel.combined.filtered.csv',
    )
    main(args.pheno_path, args.data_path, combined_path)
