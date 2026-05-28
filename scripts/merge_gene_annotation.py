#!/usr/bin/env python3
"""
merge_gene_annotation.py
------------------------
Merge gene_annotation_previous_research.csv into feature_scores_annotated.csv,
preserving all rows from feature_scores_annotated.csv (left join on feature).

Post-merge fills (each done independently):
  - gene_names_SNP1: NaN values filled from research file SNP1→gene_names lookup
  - gene_names_SNP2: NaN values filled from research file SNP2→gene_names lookup
  - previously_known: NaN values filled from per-SNP lookups; epi pairs get
    "snp1_result,snp2_result" (e.g. "yes,no")

Covariate rows ((Intercept), SEX, age, PC1–PC10) are excluded from output.
entrez_matched column is dropped.

Output: results/type2Diabetes/merged_gene_annotation.csv
"""

import argparse
import os
import re
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


def build_gene_name_lookups(df: pd.DataFrame) -> tuple[dict, dict]:
    """
    Build per-SNP gene name lookups from the research file.

    Returns
    -------
    snp1_lookup : {SNP1_value -> gene_names_SNP1}
    snp2_lookup : {SNP2_value -> gene_names_SNP2}
    """
    snp1_lookup: dict = {}
    snp2_lookup: dict = {}

    for _, row in df[df['SNP1'].notna() & df['gene_names_SNP1'].notna()].iterrows():
        key = str(row['SNP1']).strip()
        val = str(row['gene_names_SNP1']).strip()
        if key and val and val.lower() != 'nan' and key not in snp1_lookup:
            snp1_lookup[key] = val

    for _, row in df[df['SNP2'].notna() & df['gene_names_SNP2'].notna()].iterrows():
        key = str(row['SNP2']).strip()
        val = str(row['gene_names_SNP2']).strip()
        if key and val and val.lower() != 'nan' and key not in snp2_lookup:
            snp2_lookup[key] = val

    return snp1_lookup, snp2_lookup


def build_previously_known_lookup(df: pd.DataFrame) -> dict:
    """
    Build {rsid -> previously_known} from single-SNP rows in the research file.
    Used to fill missing previously_known for epi pairs and unmatched features.
    """
    single = df[df['feature'].str.match(r'^rs\d+$', case=False, na=False)]
    lookup = {}
    for _, row in single.iterrows():
        pk = str(row.get('previously_known', '')).strip()
        if pk and pk.lower() != 'nan':
            lookup[row['feature'].lower()] = pk
    return lookup


def fill_gene_names(merged: pd.DataFrame, snp1_lookup: dict, snp2_lookup: dict) -> pd.DataFrame:
    """Fill NaN gene_names_SNP1 and gene_names_SNP2 separately from per-SNP lookups."""
    merged = merged.copy()

    mask1 = merged['gene_names_SNP1'].isna()
    merged.loc[mask1, 'gene_names_SNP1'] = merged.loc[mask1, 'SNP1'].map(snp1_lookup)

    mask2 = merged['gene_names_SNP2'].isna()
    merged.loc[mask2, 'gene_names_SNP2'] = merged.loc[mask2, 'SNP2'].map(snp2_lookup)

    return merged


def fill_previously_known(merged: pd.DataFrame, snp_lookup: dict) -> pd.DataFrame:
    """
    Fill missing previously_known from per-SNP lookups.
    Epi pairs get "snp1_result,snp2_result"; single-SNP rows get the direct result.
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

    # Fill gene_names_SNP1 and gene_names_SNP2 separately
    snp1_gene_lookup, snp2_gene_lookup = build_gene_name_lookups(prev_research)
    print(f"  Gene name lookups: {len(snp1_gene_lookup)} SNP1 entries, {len(snp2_gene_lookup)} SNP2 entries")
    merged = fill_gene_names(merged, snp1_gene_lookup, snp2_gene_lookup)
    print(f"  gene_names_SNP1 filled: {merged['gene_names_SNP1'].notna().sum()} rows")
    print(f"  gene_names_SNP2 filled: {merged['gene_names_SNP2'].notna().sum()} rows")

    # Fill previously_known from per-SNP lookups
    snp_pk_lookup = build_previously_known_lookup(prev_research)
    print(f"  previously_known SNP lookup: {len(snp_pk_lookup)} entries")
    merged = fill_previously_known(merged, snp_pk_lookup)
    print(f"  {merged['previously_known'].notna().sum()} rows with previously_known after SNP fill")

    # Drop covariate rows and entrez_matched
    before = len(merged)
    merged = merged[~merged['feature'].astype(str).str.match(COVARIATE_RE)]
    print(f"\n  Dropped {before - len(merged)} covariate rows (Intercept, SEX, age, PC1-PC10)")

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
