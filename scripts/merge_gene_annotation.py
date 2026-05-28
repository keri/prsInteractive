#!/usr/bin/env python3
"""
merge_gene_annotation.py
------------------------
Merge gene_annotation_previous_research.csv into feature_scores_annotated.csv,
preserving all rows from feature_scores_annotated.csv (left join on feature).

previously_known NaN fill:
  For each rsID, check every row in the research file where that rsID appears
  as SNP1 or SNP2 (single-SNP features AND within pairs).  If the SNP is
  'yes' in any of those rows it is annotated 'yes'; otherwise 'no'.
  Epi pair rows get "snp1_result,snp2_result" (e.g. "yes,no").

Covariate rows ((Intercept), SEX, age, PC1–PC10) are excluded from output.
entrez_matched column is dropped.

Output: results/type2Diabetes/merged_gene_annotation.csv
"""

import argparse
import os
import re
from collections import defaultdict

import pandas as pd


MERGE_COLS = [
    'feature',
    'cardioFeature',
    'previously_known',
    'previous_OR_SNP1',
    'JOURNAL_SNP1',
    '95% CI (TEXT)_SNP1',
    'previous_OR_SNP2',
    'JOURNAL_SNP2',
    '95% CI (TEXT)_SNP2',
]

COVARIATE_RE = re.compile(r'^\(?intercept\)?$|^(sex|age|pc\d+)$', re.IGNORECASE)


def load_previous_research(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df = df[[c for c in df.columns if not c.startswith('Unnamed')]]
    df = df[df['feature'].notna() & (df['feature'].astype(str).str.strip() != '')]
    df['feature'] = df['feature'].astype(str).str.strip()
    return df


def deduplicate_by_feature(df: pd.DataFrame) -> pd.DataFrame:
    """One row per feature — prefer rows with previously_known set, then most data."""
    df = df.copy()
    df['_pk_filled'] = df['previously_known'].notna() & (df['previously_known'].astype(str).str.strip() != '')
    df['_n_filled']  = df.notna().sum(axis=1)
    df = (df.sort_values(['_pk_filled', '_n_filled'], ascending=[False, False])
            .drop_duplicates(subset='feature', keep='first')
            .drop(columns=['_pk_filled', '_n_filled']))
    return df


def _parse_pk(raw: str) -> tuple[str, str]:
    """
    Parse a previously_known value into (pk_snp1, pk_snp2).

    Single value ('yes' / 'no') → both positions get that value.
    Pair value ('yes,no' / 'no,yes' / 'yes,yes(BMI)') → split on first comma.
    Returns normalised 'yes' or 'no' (anything starting with 'yes' → 'yes').
    """
    raw = str(raw).strip().lower()
    if ',' in raw:
        parts = raw.split(',', 1)
        pk1 = 'yes' if parts[0].strip().startswith('yes') else ('no' if parts[0].strip().startswith('no') else '')
        pk2 = 'yes' if parts[1].strip().startswith('yes') else ('no' if parts[1].strip().startswith('no') else '')
    else:
        pk1 = 'yes' if raw.startswith('yes') else ('no' if raw.startswith('no') else '')
        pk2 = pk1
    return pk1, pk2


def build_previously_known_lookup(df: pd.DataFrame) -> dict:
    """
    Build {rsid_lower: 'yes' | 'no'} by scanning every row in the research file.

    Each rsID is checked wherever it appears as SNP1 or SNP2 — in single-SNP
    features and within pairs.  A SNP is 'yes' if it is 'yes' in any row;
    'no' only if every row that contains it says 'no'.
    """
    snp_pk_values: dict = defaultdict(set)

    for _, row in df.iterrows():
        pk_raw = str(row.get('previously_known', '')).strip()
        if not pk_raw or pk_raw.lower() == 'nan':
            continue

        pk1, pk2 = _parse_pk(pk_raw)

        snp1 = str(row.get('SNP1', '')).strip().lower()
        snp2 = str(row.get('SNP2', '')).strip().lower()

        if snp1 and snp1 != 'nan' and pk1:
            snp_pk_values[snp1].add(pk1)
        if snp2 and snp2 != 'nan' and pk2:
            snp_pk_values[snp2].add(pk2)

    # 'yes' in any row wins
    return {snp: ('yes' if 'yes' in vals else 'no')
            for snp, vals in snp_pk_values.items()}


def fill_previously_known(merged: pd.DataFrame, snp_lookup: dict) -> pd.DataFrame:
    """
    Fill missing previously_known using per-SNP lookups.
    Epi pairs: "snp1_result,snp2_result".  Single SNPs: direct result.
    """
    merged = merged.copy()
    missing_mask = (
        merged['previously_known'].isna() |
        (merged['previously_known'].astype(str).str.strip() == '') |
        (merged['previously_known'].astype(str).str.strip().str.lower() == 'nan')
    )

    for idx in merged.index[missing_mask]:
        snp1 = merged.at[idx, 'SNP1']
        snp2 = merged.at[idx, 'SNP2']

        snp1_str = str(snp1).strip().lower() if pd.notna(snp1) else ''
        snp2_str = str(snp2).strip().lower() if pd.notna(snp2) else ''

        pk1 = snp_lookup.get(snp1_str, '')
        pk2 = snp_lookup.get(snp2_str, '')

        if snp2_str and (pk1 or pk2):
            merged.at[idx, 'previously_known'] = f"{pk1},{pk2}"
        elif snp1_str and pk1:
            merged.at[idx, 'previously_known'] = pk1

    return merged


def main(pheno_path: str) -> None:
    feature_scores_path = os.path.join(pheno_path, 'feature_scores_annotated.csv')
    prev_research_path  = os.path.join(pheno_path, 'gene_annotation_previous_research.csv')
    output_path         = os.path.join(pheno_path, 'merged_gene_annotation.csv')

    print(f"Loading feature scores   : {feature_scores_path}")
    feature_scores = pd.read_csv(feature_scores_path)
    print(f"  {len(feature_scores)} rows")

    print(f"Loading previous research: {prev_research_path}")
    prev_research = load_previous_research(prev_research_path)
    print(f"  {len(prev_research)} non-blank rows, {prev_research['feature'].nunique()} unique features")

    # Merge on feature (left join)
    prev_dedup = deduplicate_by_feature(prev_research[MERGE_COLS])
    print(f"  {len(prev_dedup)} rows after deduplication")
    merged = feature_scores.merge(prev_dedup, on='feature', how='left')
    print(f"  {merged['previously_known'].notna().sum()} rows with previously_known after merge")

    # Fill missing previously_known from per-SNP lookup across all research rows
    snp_pk_lookup = build_previously_known_lookup(prev_research)
    print(f"  Per-SNP lookup built: {len(snp_pk_lookup)} unique rsIDs")
    merged = fill_previously_known(merged, snp_pk_lookup)
    print(f"  {merged['previously_known'].notna().sum()} rows with previously_known after SNP fill")

    # Drop covariate rows
    before = len(merged)
    merged = merged[~merged['feature'].astype(str).str.match(COVARIATE_RE)]
    print(f"\n  Dropped {before - len(merged)} covariate rows (Intercept, SEX, age, PC1–PC10)")

    # Drop entrez_matched
    if 'entrez_matched' in merged.columns:
        merged = merged.drop(columns=['entrez_matched'])

    merged.to_csv(output_path, index=False)
    print(f"  Wrote {len(merged)} rows to {output_path}")

    pk_counts = merged['previously_known'].value_counts(dropna=False).head(15)
    print("\npreviously_known distribution (top 15):")
    print(pk_counts.to_string())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge gene annotation with previous research")
    parser.add_argument(
        '--pheno_path',
        default='/Users/kerimulterer/prsInteractive/results/type2Diabetes',
        help='Path to the phenotype results directory',
    )
    args = parser.parse_args()
    main(args.pheno_path)
