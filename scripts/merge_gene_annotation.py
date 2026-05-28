#!/usr/bin/env python3
"""
merge_gene_annotation.py
------------------------
Merge gene_annotation_previous_research.csv into feature_scores_annotated.csv,
preserving all rows from feature_scores_annotated.csv (left join on feature).

Missing previously_known values are filled from individual SNP lookups:
  - For epi pairs (SNP1,SNP2): look up each SNP separately and combine as
    "snp1_result,snp2_result" (e.g. "yes,no").
  - For single-SNP features: use the SNP's own lookup result directly.

Output: results/type2Diabetes/merged_gene_annotation.csv
"""

import argparse
import os
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


def load_previous_research(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    # Drop Excel padding columns and blank rows
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


def build_snp_lookup(df: pd.DataFrame) -> dict:
    """
    Build {rsid_lower: previously_known} from single-SNP rows in the research file.
    Used to fill missing previously_known for epi pairs.
    """
    single = df[df['feature'].str.match(r'^rs\d+$', case=False, na=False)]
    lookup = {}
    for _, row in single.iterrows():
        pk = str(row.get('previously_known', '')).strip()
        if pk and pk.lower() != 'nan':
            lookup[row['feature'].lower()] = pk
    return lookup


def fill_previously_known(merged: pd.DataFrame, snp_lookup: dict) -> pd.DataFrame:
    """
    For rows with missing previously_known, look up SNP1 and SNP2 individually
    and combine as "snp1_result,snp2_result".
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
            # epi pair — combine as "snp1_result,snp2_result"
            merged.at[idx, 'previously_known'] = f"{pk1},{pk2}"
        elif snp1_str and pk1:
            # single SNP
            merged.at[idx, 'previously_known'] = pk1

    return merged


def main(pheno_path: str) -> None:
    feature_scores_path   = os.path.join(pheno_path, 'feature_scores_annotated.csv')
    prev_research_path    = os.path.join(pheno_path, 'gene_annotation_previous_research.csv')
    output_path           = os.path.join(pheno_path, 'merged_gene_annotation.csv')

    print(f"Loading feature scores  : {feature_scores_path}")
    feature_scores = pd.read_csv(feature_scores_path)
    print(f"  {len(feature_scores)} rows")

    print(f"Loading previous research: {prev_research_path}")
    prev_research = load_previous_research(prev_research_path)
    print(f"  {len(prev_research)} non-blank rows, {prev_research['feature'].nunique()} unique features")

    prev_dedup = deduplicate_by_feature(prev_research[MERGE_COLS])
    print(f"  {len(prev_dedup)} rows after deduplication")

    # Task 2: left merge preserving all feature_scores rows
    merged = feature_scores.merge(prev_dedup, on='feature', how='left')
    print(f"  {merged['previously_known'].notna().sum()} rows with previously_known after merge")

    # Task 3: fill missing previously_known from individual SNP lookups
    snp_lookup = build_snp_lookup(prev_research)
    print(f"  SNP lookup built: {len(snp_lookup)} entries")
    merged = fill_previously_known(merged, snp_lookup)
    filled = merged['previously_known'].notna().sum()
    print(f"  {filled} rows with previously_known after SNP fill")

    # Task 4: save
    merged.to_csv(output_path, index=False)
    print(f"\nWrote {len(merged)} rows to {output_path}")

    # Summary of previously_known distribution
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
