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
  high_risk_cases: PHENOTYPE == 2  AND  top-two-decile PRS bin (bin == 9 or 10)

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
  <pheno_path>/env_summary_by_prs_group.csv               — all E features
  <pheno_path>/env_summary_by_prs_group.holdout.csv       — all E features (holdout)
  <pheno_path>/env_summary_cardio_features.csv            — cardio E terms only (table-ready)
  <pheno_path>/env_summary_cardio_features.holdout.csv    — cardio E terms only (holdout)
  <pheno_path>/env_summary_clinical_measures.holdout.csv  — clinical measure means + binary
                                                            prevalence from ClinicalMeasures file
                                                            (PHENOTYPE 0/1; high-risk via
                                                            {model}_high_risk flag)
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
    'bin_combined':       'combined',
}

SCALED_TO_BIN = {
    'scaled_prs_epi_product':    'bin_epi_product',
    'scaled_prs_cardio_product': 'bin_cardio_product',
    'scaled_prs_main':           'bin_main',
    'scaled_prs_epi_summed':     'bin_epi_summed',
    'scaled_prs_cardio_summed':  'bin_cardio_summed',
    'scaled_prs_combined':       'bin_combined',
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


def summarise(df: pd.DataFrame, env_cols: list, diagnosis_code: str = '') -> pd.DataFrame:
    """
    diagnosis_code: 'E11' for T2D, 'K90' for celiac.  Only that code's
    pct_after_assessment column is written; pass '' to compute both.
    The pct is per-person: fraction of individuals (among those with a
    known diagnosis date) whose age-at-diagnosis > age-at-participation.
    """
    records = []

    for bin_col, model in BIN_TO_MODEL.items():
        if bin_col not in df.columns:
            continue

        # Top-two-decile threshold: works for both 1-10 and 1-1000 bin scales
        bin_max = df[bin_col].max()
        high_risk_threshold = bin_max * 0.8

        groups = {
            'controls':        df[df['PHENOTYPE'] == 1],
            'all_cases':       df[df['PHENOTYPE'] == 2],
            'high_risk_cases': df[(df['PHENOTYPE'] == 2) & (df[bin_col] > high_risk_threshold)],
        }

        codes_to_report = (diagnosis_code,) if diagnosis_code else ('E11', 'K90')

        for group_name, grp in groups.items():
            rec: dict = {'model': model, 'group': group_name, 'n': len(grp)}

            rec['pct_male']                  = _pct_male(grp.get('is_male', pd.Series(dtype=float)))
            rec['age_at_participation_mean'] = _mean(grp.get('age_at_participation', pd.Series(dtype=float)))
            rec['age_E11_diagnosis_mean']    = _mean(grp.get('age_E11_diagnosis', pd.Series(dtype=float)))
            rec['age_K90_diagnosis_mean']    = _mean(grp.get('age_K90_diagnosis', pd.Series(dtype=float)))

            # Per-person: fraction of individuals (with known diagnosis date)
            # whose age-at-diagnosis > age-at-recruitment.
            part_col = 'age_at_participation'
            for code in codes_to_report:
                diag_col = f'age_{code}_diagnosis'
                if diag_col in grp.columns and part_col in grp.columns:
                    both_known = grp[[diag_col, part_col]].dropna()
                    if len(both_known) > 0:
                        pct = round(
                            (both_known[diag_col] > both_known[part_col]).mean() * 100, 1
                        )
                        rec[f'pct_{code}_diagnosis_after_assessment'] = pct

            for feat in env_cols:
                if feat in grp.columns:
                    rec[f'{feat}_mean'] = _mean(grp[feat])

            records.append(rec)

    return pd.DataFrame(records)


def process_file(
    groups_path: str,
    demog: pd.DataFrame,
    env: pd.DataFrame,
    diagnosis_code: str = '',
) -> pd.DataFrame:
    df = pd.read_csv(groups_path)
    df = add_bins_from_scaled(df)

    env_renamed = env.rename(columns={'Participant ID': 'IID'})
    df = df.merge(demog,        on='IID', how='left')
    df = df.merge(env_renamed,  on='IID', how='left')

    env_cols = [c for c in env.columns if c != 'Participant ID']
    return summarise(df, env_cols, diagnosis_code=diagnosis_code)


SKIP_PREFIXES = (
    'scaled_prs', 'prs_', 'bin_', 'any_', 'n_', 'prop_',
    'combined', 'PHENOTYPE', 'IID',
)
SKIP_SUFFIXES = ('_case_bin', '_high_risk', '_threshold_value', '_binary')


def load_clinical_binary_flags(path: str) -> tuple[list[str], pd.DataFrame]:
    """
    Load the ClinicalMeasures file.

    Returns
    -------
    features : list[str]
        Continuous feature names (= binary col names with '_binary' stripped).
        These match column names in participant_environment.csv.
    df : pd.DataFrame
        IID, PHENOTYPE (remapped 0→1 / 1→2), {model}_high_risk flags,
        {feature}_binary columns.
    """
    df = pd.read_csv(path)
    df['PHENOTYPE'] = df['PHENOTYPE'].map({0: 1, 1: 2})

    binary_cols = [c for c in df.columns if c.endswith('_binary')]
    hr_cols     = [c for c in df.columns if c.endswith('_high_risk')]
    features    = [c.replace('_binary', '') for c in binary_cols]

    keep = ['IID', 'PHENOTYPE'] + hr_cols + binary_cols
    return features, df[keep]


def summarise_clinical_binary(
    clinical_df: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    """
    Compute binary prevalence (% positive) by model / group from the
    ClinicalMeasures file.  High-risk defined by {model}_high_risk == 1
    (this threshold differs from the top-decile bin==10 used in the
    main env summary).
    """
    records = []
    hr_cols = [c for c in clinical_df.columns if c.endswith('_high_risk')]

    for hr_col in hr_cols:
        model = hr_col.replace('_high_risk', '')
        groups = {
            'controls':        clinical_df[clinical_df['PHENOTYPE'] == 1],
            'all_cases':       clinical_df[clinical_df['PHENOTYPE'] == 2],
            'high_risk_cases': clinical_df[
                (clinical_df['PHENOTYPE'] == 2) & (clinical_df[hr_col] == 1)
            ],
        }
        for group_name, grp in groups.items():
            rec: dict = {'model': model, 'group': group_name,
                         'n_clinical': len(grp)}
            for feat in features:
                bin_col = f'{feat}_binary'
                if bin_col in grp.columns:
                    v = grp[bin_col].dropna()
                    rec[f'{feat}_pct_binary'] = (
                        round(float(v.mean() * 100), 1) if len(v) > 0 else None
                    )
            records.append(rec)

    return pd.DataFrame(records)


def make_clinical_comparison(
    disc_summary: pd.DataFrame,
    hold_summary: pd.DataFrame,
    hold_binary: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    """
    Build a side-by-side discovery / holdout comparison for clinical E features.

    Continuous means come from the two env summaries (participant_environment.csv
    source, consistent bin==10 high-risk definition).
    Binary prevalence comes from hold_binary (ClinicalMeasures file, uses the
    {model}_high_risk flag which may cover a broader high-risk set).

    Returns a wide DataFrame with one row per (model, group).
    """
    mean_cols = {f: f'{f}_mean' for f in features}

    def _extract(summary: pd.DataFrame, prefix: str) -> pd.DataFrame:
        cols = ['model', 'group'] + [
            c for c in summary.columns if c in mean_cols.values()
        ]
        sub = summary[cols].copy()
        sub = sub.rename(columns={v: f'{prefix}_{v}' for v in mean_cols.values()})
        return sub

    disc = _extract(disc_summary, 'discovery')
    hold = _extract(hold_summary, 'holdout')

    merged = disc.merge(hold, on=['model', 'group'], how='outer')
    merged = merged.merge(hold_binary, on=['model', 'group'], how='left')
    return merged


def get_cardio_e_features(combined_path: str, env_cols: set[str]) -> list[str]:
    """Return sorted list of E feature names with non-zero betas in cardio
    (GxGxE / G+GxE) models that are present in participant_environment.csv."""
    combined = pd.read_csv(combined_path)
    cardio = combined[
        combined['model'].isin({'cardio_product', 'cardio_summed'}) &
        (combined['coefs'] != 0)
    ]
    features: set[str] = set()
    for f in cardio['feature'].astype(str):
        if ',' in f:
            tok = f.split(',')[0]
            if tok not in SKIP_TERMS and not _RS_RE.match(tok) and tok in env_cols:
                features.add(tok)
    return sorted(features)


def main(pheno_path: str, data_path: str, combined_path: str) -> None:
    # Detect which diagnosis code is the primary outcome for this phenotype
    path_lower = pheno_path.lower()
    if 'type2' in path_lower or 't2d' in path_lower:
        diagnosis_code = 'E11'
    elif 'celiac' in path_lower or 'coeliac' in path_lower:
        diagnosis_code = 'K90'
    else:
        diagnosis_code = ''
    print(f"Primary diagnosis code: {diagnosis_code or '(both E11 and K90)'}")

    print(f"Loading demographics : {data_path}/participant.csv")
    demog = load_demographics(data_path)
    print(f"  {len(demog)} participants")

    env_path = os.path.join(data_path, 'participant_environment.csv')
    print(f"Loading environment  : {env_path}")
    env = pd.read_csv(env_path)
    print(f"  {len(env)} participants, {len(env.columns) - 1} E features")

    env_cols = set(env.columns) - {'Participant ID'}
    print(f"Loading cardio E terms: {combined_path}")
    cardio_features = get_cardio_e_features(combined_path, env_cols)
    print(f"  {len(cardio_features)} unique E features (non-zero beta, cardio models): {cardio_features}")

    # Load clinical measures file (holdout only)
    clinical_path = os.path.join(
        pheno_path, 'combinedAnalysis', 'scores',
        'clinicalMeasures', 'combinedPRS.holdout.ClinicalMeasures.csv',
    )
    clinical_features: list[str] = []
    hold_binary_summary: pd.DataFrame = pd.DataFrame()
    if os.path.exists(clinical_path):
        print(f"Loading clinical measures: {clinical_path}")
        clinical_features, clinical_flags = load_clinical_binary_flags(clinical_path)
        print(f"  Clinical E features: {clinical_features}")
        hold_binary_summary = summarise_clinical_binary(clinical_flags, clinical_features)
    else:
        print(f"[skip] clinical measures not found: {clinical_path}")

    scores_dir = os.path.join(pheno_path, 'combinedAnalysis', 'scores')
    summaries: dict[str, pd.DataFrame] = {}

    for fname, label, out_name, cardio_out_name in [
        ('combinedPRSGroups.filtered.csv',         'discovery',
         'env_summary_by_prs_group.csv',
         'env_summary_cardio_features.csv'),
        ('combinedPRSGroups.holdout.filtered.csv', 'holdout',
         'env_summary_by_prs_group.holdout.csv',
         'env_summary_cardio_features.holdout.csv'),
    ]:
        fpath = os.path.join(scores_dir, fname)
        if not os.path.exists(fpath):
            print(f"\n[skip] not found: {fpath}")
            continue

        print(f"\nProcessing {label}: {fpath}")
        summary = process_file(fpath, demog, env, diagnosis_code=diagnosis_code)
        summaries[label] = summary

        # Full E-feature output
        out_path = os.path.join(pheno_path, out_name)
        summary.to_csv(out_path, index=False)
        print(f"  Wrote {len(summary)} rows (all E features) → {out_path}")

        # Cardio-features-only table output — include only the primary diagnosis pct column
        diag_pct_col = f'pct_{diagnosis_code}_diagnosis_after_assessment' if diagnosis_code else None
        meta_cols = ['model', 'group', 'n', 'pct_male',
                     'age_at_participation_mean',
                     'age_E11_diagnosis_mean', 'age_K90_diagnosis_mean']
        if diag_pct_col:
            # Insert pct column immediately after the relevant diagnosis age column
            age_col = f'age_{diagnosis_code}_diagnosis_mean'
            idx = meta_cols.index(age_col) + 1 if age_col in meta_cols else len(meta_cols)
            meta_cols.insert(idx, diag_pct_col)
        else:
            meta_cols += ['pct_E11_diagnosis_after_assessment',
                          'pct_K90_diagnosis_after_assessment']
        cardio_mean_cols = [f'{f}_mean' for f in cardio_features if f'{f}_mean' in summary.columns]
        keep_cols = [c for c in meta_cols if c in summary.columns] + cardio_mean_cols
        cardio_out_path = os.path.join(pheno_path, cardio_out_name)
        summary[keep_cols].to_csv(cardio_out_path, index=False)
        print(f"  Wrote {len(summary)} rows ({len(cardio_mean_cols)} cardio E features) → {cardio_out_path}")

        preview_cols = ['model', 'group', 'n', 'pct_male',
                        'age_at_participation_mean',
                        'age_E11_diagnosis_mean', 'age_K90_diagnosis_mean']
        preview_cols = [c for c in preview_cols if c in summary.columns]
        print(summary[preview_cols].to_string(index=False))

    # Clinical measures comparison: discovery vs holdout continuous means + holdout binary %
    if clinical_features and 'discovery' in summaries and 'holdout' in summaries:
        comparison = make_clinical_comparison(
            summaries['discovery'], summaries['holdout'],
            hold_binary_summary, clinical_features,
        )
        clin_out = os.path.join(pheno_path, 'env_summary_clinical_measures.csv')
        comparison.to_csv(clin_out, index=False)
        print(f"\nClinical measures comparison → {clin_out}")
        # Preview: show mean cols and binary % for high_risk_cases
        preview = comparison[comparison['group'] == 'high_risk_cases'].copy()
        mean_c = [c for c in comparison.columns if c.endswith('_mean')]
        bin_c  = [c for c in comparison.columns if c.endswith('_pct_binary')]
        print(preview[['model', 'group', 'n_clinical'] + mean_c + bin_c].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarise E features by PRS risk group for coauthor review"
    )
    parser.add_argument(
        '--pheno_path', required=True,
        default='/Users/kerimulterer/prsInteractive/results/type2Diabetes',
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