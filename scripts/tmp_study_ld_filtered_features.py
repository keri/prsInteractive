#!/usr/bin/env python3
"""
Temporary study script: compare importantFeaturesPostShap.filteredLD.csv
(output of run_plink_LD.sh) against the pre-LD input for both epi subfolders.

Run after:
  ./run_plink_LD.sh type2Diabetes productEpi
  ./run_plink_LD.sh type2Diabetes summedEpi

Output: results/type2Diabetes/ldFilterStudy.csv
"""

import pandas as pd
import os

PHENO_PATH = "/Users/kerimulterer/prsInteractive/results/type2Diabetes"

subfolders = ["productEpi", "summedEpi"]
records = []

for sub in subfolders:
    scores_dir = f"{PHENO_PATH}/{sub}/scores"
    pre_ld  = f"{scores_dir}/importantFeaturesPostShap.csv"
    post_ld = f"{scores_dir}/importantFeaturesPostShap.filteredLD.csv"

    if not os.path.exists(post_ld):
        print(f"[SKIP] {sub}: filteredLD file not found — run run_plink_LD.sh {sub} first")
        continue

    df_pre  = pd.read_csv(pre_ld)
    df_post = pd.read_csv(post_ld)

    pre_features  = set(df_pre['feature'])
    post_features = set(df_post['feature'])
    pruned        = pre_features - post_features

    # tag each post-LD row with source subfolder
    df_post = df_post.copy()
    df_post['subfolder'] = sub
    df_post['ld_pruned'] = False
    records.append(df_post)

    # also add pruned rows so we can see what was removed
    df_pruned = df_pre[df_pre['feature'].isin(pruned)].copy()
    df_pruned['subfolder'] = sub
    df_pruned['ld_pruned'] = True
    records.append(df_pruned)

    print(f"\n{sub}")
    print(f"  pre-LD features : {len(pre_features)}")
    print(f"  post-LD features: {len(post_features)}")
    print(f"  pruned by LD    : {len(pruned)}")
    for dt in df_pre['data_type'].unique():
        n_pre    = (df_pre['data_type']  == dt).sum()
        n_post   = (df_post[df_post['ld_pruned'] == False]['data_type'] == dt).sum()
        n_pruned = n_pre - n_post
        print(f"    {dt:12s}  {n_pre:4d} -> {n_post:4d}  ({n_pruned} pruned)")

if not records:
    print("\nNo filteredLD files found. Run run_plink_LD.sh for each subfolder first.")
else:
    out = pd.concat(records, ignore_index=True)
    out_path = f"{PHENO_PATH}/ldFilterStudy.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} rows to {out_path}")
