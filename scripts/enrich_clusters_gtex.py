#!/usr/bin/env python3
"""
enrich_clusters_gtex.py
=======================
Post-NMF enrichment: query GTEx for tissue-specific expression and eQTL evidence
for the top-loaded genomic features in each cluster.

Workflow
--------
1. Load feature_loadings.csv from a bNMF output directory.
2. Extract top N gen_* features per cluster.
3. Strip prefix → SNP rs-IDs or interaction pairs.
4. Map SNPs → gene symbols via NCBI dbSNP E-utilities.
5. Query GTEx REST API:
   - Tissue-specific median expression per gene.
   - eQTL effect sizes per SNP × gene × tissue.
6. Score: gtex_score = |eqtl_effect_size| × median_TPM
7. Write enrichment tables and heatmap to {cluster_dir}/enrichment/gtex/

Usage
-----
  python enrich_clusters_gtex.py \
    --cluster_dir /path/to/combinedCohortBnmf \
    --phenotypes "type 2 diabetes,epilepsy" \
    --top_n 20 \
    --ncbi_api_key YOUR_KEY
"""

import argparse
import json
import os
import time
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns


# ── Constants ─────────────────────────────────────────────────────────────────

GTEX_API    = "https://gtexportal.org/rest/v1"
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GENOME_BUILD = "GRCh38"

# Tissues most relevant to cardiometabolic / neurological phenotypes.
# Used as the default subset when --tissues is not specified.
DEFAULT_TISSUES = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Hippocampus",
    "Heart_Left_Ventricle", "Heart_Atrial_Appendage",
    "Kidney_Cortex",
    "Liver",
    "Lung",
    "Muscle_Skeletal",
    "Pancreas",
    "Small_Intestine_Terminal_Ileum",
    "Whole_Blood",
]

CACHE_DIR   = None   # set at runtime


# ── HTTP helpers ───────────────────────────────────────────────────────────────

def _get(url, params=None, retries=3, wait=1.5, headers=None):
    """GET with retry on rate-limit or transient errors."""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=30)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429:
                time.sleep(wait * (attempt + 1))
                continue
            r.raise_for_status()
        except requests.RequestException as e:
            if attempt == retries - 1:
                print(f"    [WARN] Request failed after {retries} retries: {e}")
                return None
            time.sleep(wait)
    return None


def _ncbi_headers(api_key):
    h = {'User-Agent': 'enrich_clusters_gtex/1.0 (research; contact via GitHub)'}
    if api_key:
        h['api_key'] = api_key
    return h


# ── SNP → gene mapping ────────────────────────────────────────────────────────

def _cache_path(name):
    return os.path.join(CACHE_DIR, f"{name}.json") if CACHE_DIR else None


def _load_cache(name):
    p = _cache_path(name)
    if p and os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return {}


def _save_cache(name, data):
    p = _cache_path(name)
    if p:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'w') as f:
            json.dump(data, f)


def snps_to_genes(snp_ids, ncbi_api_key=None):
    """
    Map a list of rs-IDs to gene symbols via NCBI E-utilities.

    Returns
    -------
    dict  {rs_id: [gene_symbol, ...]}
    """
    cache = _load_cache('snp_gene_map')
    to_fetch = [s for s in snp_ids if s not in cache and s.startswith('rs')]
    headers  = _ncbi_headers(ncbi_api_key)

    batch_size = 200
    rate_limit = 0.12 if ncbi_api_key else 0.35   # seconds between requests

    for i in range(0, len(to_fetch), batch_size):
        batch = to_fetch[i:i + batch_size]
        ids   = ','.join(b.lstrip('rs') for b in batch)
        data  = _get(
            f"{NCBI_BASE}/esummary.fcgi",
            params={'db': 'snp', 'id': ids, 'retmode': 'json'},
            headers=headers,
        )
        if data and 'result' in data:
            for rs_num, rec in data['result'].items():
                if rs_num == 'uids':
                    continue
                genes = [g.get('name', '') for g in rec.get('genes', [])]
                genes = [g for g in genes if g]
                cache[f'rs{rs_num}'] = genes
        time.sleep(rate_limit)

    _save_cache('snp_gene_map', cache)
    return {s: cache.get(s, []) for s in snp_ids}


def parse_feature_to_snps(feature_name):
    """
    Strip gen_ prefix and split interaction pairs into component SNP rs-IDs.
    Returns list of rs-IDs (may be empty for non-rs features).
    """
    name = feature_name.removeprefix('gen_')
    parts = name.split('_') if '_' in name else [name]
    return [p for p in parts if p.startswith('rs')]


# ── GTEx queries ──────────────────────────────────────────────────────────────

def get_gencode_id(gene_symbol):
    """Return GENCODE ID (ENSG...) for a gene symbol via GTEx API."""
    data = _get(
        f"{GTEX_API}/reference/gene",
        params={'geneName': gene_symbol, 'referenceGenomeBuild': GENOME_BUILD},
    )
    if data and data.get('gene'):
        return data['gene'][0].get('gencodeId')
    return None


def get_tissue_expression(gencode_id, tissues):
    """
    Fetch median TPM per tissue for one gene.

    Returns
    -------
    dict  {tissue_id: median_tpm}
    """
    results = {}
    for tissue in tissues:
        data = _get(
            f"{GTEX_API}/expression/geneExpression",
            params={
                'gencodeId':         gencode_id,
                'tissueSiteDetailId': tissue,
            },
        )
        if data and data.get('geneExpression'):
            entry  = data['geneExpression'][0]
            tpm    = entry.get('median', entry.get('data', [None])[0])
            if tpm is not None:
                results[tissue] = float(tpm)
        time.sleep(0.2)
    return results


def get_eqtl_evidence(snp_id, gencode_id, tissues):
    """
    Fetch eQTL effect sizes for one SNP × gene pair across tissues.

    Returns
    -------
    dict  {tissue_id: {'effect_size': float, 'p_value': float}}
    """
    results = {}
    for tissue in tissues:
        data = _get(
            f"{GTEX_API}/eqtl/singleTissueEqtl",
            params={
                'variantId':          snp_id,
                'gencodeId':          gencode_id,
                'tissueSiteDetailId': tissue,
            },
        )
        if data and data.get('singleTissueEqtl'):
            rec = data['singleTissueEqtl'][0]
            results[tissue] = {
                'effect_size': float(rec.get('nes', 0)),
                'p_value':     float(rec.get('pValue', 1.0)),
            }
        time.sleep(0.2)
    return results


# ── Main enrichment ───────────────────────────────────────────────────────────

def enrich_with_gtex(
    cluster_dir,
    top_n=20,
    tissues=None,
    ncbi_api_key=None,
):
    """
    Main enrichment routine.

    Parameters
    ----------
    cluster_dir   : str — bNMF output directory containing feature_loadings.csv
    top_n         : int — top features per cluster to query
    tissues       : list[str] | None — GTEx tissue IDs; None = DEFAULT_TISSUES
    ncbi_api_key  : str | None

    Returns
    -------
    dict  {'expression_df', 'eqtl_df', 'cluster_tissue_summary'}
    """
    global CACHE_DIR
    CACHE_DIR = os.path.join(cluster_dir, 'enrichment', '.cache')

    out_dir = os.path.join(cluster_dir, 'enrichment', 'gtex')
    os.makedirs(out_dir, exist_ok=True)

    if tissues is None:
        tissues = DEFAULT_TISSUES

    # ── Load feature loadings ──────────────────────────────────────────────
    loadings_path = os.path.join(cluster_dir, 'feature_loadings.csv')
    if not os.path.exists(loadings_path):
        print(f"  [ERROR] feature_loadings.csv not found in {cluster_dir}")
        return {}

    H_df = pd.read_csv(loadings_path, index_col=0)
    print(f"  Loaded feature loadings: {H_df.shape[0]} clusters × {H_df.shape[1]} features")

    # ── Collect top gen_ features per cluster ─────────────────────────────
    gen_cols = [c for c in H_df.columns if 'gen_' in c]
    all_top_features: set[str] = set()
    cluster_top: dict[str, list[str]] = {}

    for cluster in H_df.index:
        row  = H_df.loc[cluster, gen_cols]
        top  = row.nlargest(top_n).index.tolist()
        cluster_top[cluster] = top
        all_top_features.update(top)

    print(f"  Unique genomic features to query: {len(all_top_features)}")

    # ── Map SNPs → genes ───────────────────────────────────────────────────
    snp_ids: set[str] = set()
    feat_to_snps: dict[str, list[str]] = {}
    for feat in all_top_features:
        snps = parse_feature_to_snps(feat)
        feat_to_snps[feat] = snps
        snp_ids.update(snps)

    print(f"  Mapping {len(snp_ids)} SNPs to genes via NCBI...")
    snp_gene_map = snps_to_genes(list(snp_ids), ncbi_api_key=ncbi_api_key)

    # All unique genes
    all_genes: set[str] = set()
    for genes in snp_gene_map.values():
        all_genes.update(genes)
    print(f"  Genes found: {len(all_genes)}")

    # ── Fetch GTEx expression ──────────────────────────────────────────────
    print(f"\n  Fetching tissue expression for {len(all_genes)} genes "
          f"across {len(tissues)} tissues...")

    expr_cache = _load_cache('gtex_expression')
    expr_rows  = []

    for gene in sorted(all_genes):
        if gene in expr_cache:
            tpm_map = expr_cache[gene]
        else:
            gencode = get_gencode_id(gene)
            if gencode is None:
                print(f"    [SKIP] {gene}: not found in GTEx reference")
                continue
            tpm_map = get_tissue_expression(gencode, tissues)
            expr_cache[gene] = tpm_map
            _save_cache('gtex_expression', expr_cache)

        for tissue, tpm in tpm_map.items():
            expr_rows.append({'gene': gene, 'tissue': tissue, 'median_tpm': tpm})

    expr_df = pd.DataFrame(expr_rows)
    if not expr_df.empty:
        expr_df.to_csv(os.path.join(out_dir, 'gtex_gene_expression.csv'), index=False)
        print(f"    Saved gtex_gene_expression.csv ({len(expr_df)} rows)")

    # ── Fetch eQTL evidence ────────────────────────────────────────────────
    print(f"\n  Fetching eQTL evidence for {len(snp_ids)} SNPs...")

    eqtl_cache = _load_cache('gtex_eqtl')
    eqtl_rows  = []

    for snp in sorted(snp_ids):
        genes = snp_gene_map.get(snp, [])
        for gene in genes:
            cache_key = f'{snp}__{gene}'
            if cache_key in eqtl_cache:
                eqtl_map = eqtl_cache[cache_key]
            else:
                gencode = get_gencode_id(gene)
                if gencode is None:
                    continue
                eqtl_map = get_eqtl_evidence(snp, gencode, tissues)
                eqtl_cache[cache_key] = eqtl_map
                _save_cache('gtex_eqtl', eqtl_cache)

            for tissue, ev in eqtl_map.items():
                eqtl_rows.append({
                    'snp':         snp,
                    'gene':        gene,
                    'tissue':      tissue,
                    'effect_size': ev['effect_size'],
                    'p_value':     ev['p_value'],
                })

    eqtl_df = pd.DataFrame(eqtl_rows)
    if not eqtl_df.empty:
        eqtl_df.to_csv(os.path.join(out_dir, 'gtex_eqtl_evidence.csv'), index=False)
        print(f"    Saved gtex_eqtl_evidence.csv ({len(eqtl_df)} rows)")

    # ── Cluster × tissue summary ───────────────────────────────────────────
    summary_rows = []

    for cluster, top_feats in cluster_top.items():
        tissue_scores: defaultdict = defaultdict(list)

        for feat in top_feats:
            loading = float(H_df.loc[cluster, feat])
            snps    = feat_to_snps.get(feat, [])
            for snp in snps:
                genes = snp_gene_map.get(snp, [])
                for gene in genes:
                    # eQTL-weighted expression score
                    eqtl_for_snp_gene = eqtl_df[
                        (eqtl_df['snp'] == snp) & (eqtl_df['gene'] == gene)
                    ] if not eqtl_df.empty else pd.DataFrame()

                    gene_expr = expr_df[expr_df['gene'] == gene] \
                        if not expr_df.empty else pd.DataFrame()

                    for tissue in tissues:
                        tpm = gene_expr.loc[gene_expr['tissue'] == tissue, 'median_tpm']
                        eff = eqtl_for_snp_gene.loc[
                            eqtl_for_snp_gene['tissue'] == tissue, 'effect_size'
                        ]
                        tpm_val = float(tpm.iloc[0]) if not tpm.empty else 0.0
                        eff_val = abs(float(eff.iloc[0])) if not eff.empty else 0.0
                        score   = loading * eff_val * tpm_val
                        if score > 0:
                            tissue_scores[tissue].append(score)

        for tissue, scores in tissue_scores.items():
            summary_rows.append({
                'cluster':         cluster,
                'tissue':          tissue,
                'mean_gtex_score': round(np.mean(scores), 6),
                'n_genes':         len(scores),
            })

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df.sort_values(['cluster', 'mean_gtex_score'],
                               ascending=[True, False], inplace=True)
        summary_df.to_csv(
            os.path.join(out_dir, 'gtex_cluster_tissue_summary.csv'), index=False
        )
        print(f"    Saved gtex_cluster_tissue_summary.csv")
        _plot_cluster_tissue_heatmap(summary_df, out_dir)

    print(f"\n  GTEx enrichment complete → {out_dir}/")
    return {
        'expression_df':        expr_df,
        'eqtl_df':              eqtl_df,
        'cluster_tissue_summary': summary_df,
    }


def _plot_cluster_tissue_heatmap(summary_df, out_dir):
    if summary_df.empty:
        return

    pivot = summary_df.pivot_table(
        index='tissue', columns='cluster',
        values='mean_gtex_score', aggfunc='mean',
    ).fillna(0)

    # Keep top 20 tissues by max score
    if len(pivot) > 20:
        top_tissues = pivot.max(axis=1).nlargest(20).index
        pivot = pivot.loc[top_tissues]

    fig_h = max(5, len(pivot) * 0.45)
    fig_w = max(4, len(pivot.columns) * 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        pivot,
        cmap='YlOrRd',
        annot=True, fmt='.3f', annot_kws={'size': 7},
        linewidths=0.3,
        ax=ax,
        cbar_kws={'label': 'Mean GTEx score\n(loading × |eQTL| × median TPM)'},
    )
    ax.set_title('Cluster × Tissue GTEx Enrichment', fontsize=10)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('GTEx Tissue')
    plt.tight_layout()

    path = os.path.join(out_dir, 'gtex_cluster_tissue_heatmap.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved gtex_cluster_tissue_heatmap.png")


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Post-NMF GTEx enrichment: map top cluster features to tissue expression "
            "and eQTL evidence via GTEx REST API + NCBI E-utilities."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '--cluster_dir', required=True,
        help="Path to bNMF output directory containing feature_loadings.csv",
    )
    parser.add_argument(
        '--top_n', type=int, default=20,
        help="Top features per cluster to query (default: 20)",
    )
    parser.add_argument(
        '--tissues', default=None,
        help="Comma-separated GTEx tissue IDs (default: cardiometabolic panel)",
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help="NCBI Entrez API key for higher rate limits (10/sec vs 3/sec)",
    )

    args = parser.parse_args()

    tissues = [t.strip() for t in args.tissues.split(',')] \
        if args.tissues else None

    enrich_with_gtex(
        cluster_dir=args.cluster_dir,
        top_n=args.top_n,
        tissues=tissues,
        ncbi_api_key=args.ncbi_api_key,
    )
