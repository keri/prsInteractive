#!/usr/bin/env python3
"""
prepare_snpToGene_input.py
--------------------------
Merge importantFeaturesPostShap.filteredLD.csv from productEpi/ and summedEpi/,
remap data_type → PRS_calculation, deduplicate G (main) features.

PRS_calculation mapping:
  productEpi: epi → GxG, main → G
  summedEpi:  epi → G+G, main → G

G features are identical in both subfolders; productEpi copy is kept.

Output: results/type2Diabetes/snpToGeneInput.csv
"""

import argparse
import os

import pandas as pd

PRODUCT_EPI_MAP = {'epi': 'GxG', 'main': 'G'}
SUMMED_EPI_MAP  = {'epi': 'G+G', 'main': 'G'}


def load_filtered_ld(pheno_path: str, subfolder: str, dtype_map: dict) -> pd.DataFrame:
    path = os.path.join(pheno_path, subfolder, 'scores', 'importantFeaturesPostShap.filteredLD.csv')
    print(f"Loading: {path}")
    df = pd.read_csv(path)
    df['PRS_calculation'] = df['data_type'].map(dtype_map)
    df = df.drop(columns=[c for c in ['data_type', 'iteration'] if c in df.columns])
    df = df[['feature', 'shap_zscore', 'PRS_calculation']]
    print(f"  {len(df)} rows: {df['PRS_calculation'].value_counts().to_dict()}")
    return df


def main(pheno_path: str) -> None:
    product = load_filtered_ld(pheno_path, 'productEpi', PRODUCT_EPI_MAP)
    summed  = load_filtered_ld(pheno_path, 'summedEpi',  SUMMED_EPI_MAP)

    g_features     = product[product['PRS_calculation'] == 'G']
    gxg_features   = product[product['PRS_calculation'] == 'GxG']
    gplus_features = summed[summed['PRS_calculation'] == 'G+G']

    combined = pd.concat([g_features, gxg_features, gplus_features], ignore_index=True)
    print(f"\nCombined: {len(combined)} rows")
    print(combined['PRS_calculation'].value_counts().to_string())

    out_path = os.path.join(pheno_path, 'snpToGeneInput.csv')
    combined.to_csv(out_path, index=False)
    print(f"\nWrote {len(combined)} rows to {out_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare snpToGene.py input from filteredLD files")
    parser.add_argument(
        '--pheno_path',
        default='/Users/kerimulterer/prsInteractive/results/type2Diabetes',
        help='Path to phenotype results directory',
    )
    args = parser.parse_args()
    main(args.pheno_path)
