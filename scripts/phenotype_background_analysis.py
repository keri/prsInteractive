#!/usr/bin/env python3
"""
phenotype_background_analysis.py
=================================
Pre-enrichment background knowledge table for a phenotype of interest.

Uses a manually curated seed file (phenotype_seeds/<phenotype>.json) that
lists known subtypes, co-phenotypes, molecular mechanisms, pathways,
complications, comorbidities and protective factors.  Each seed term is
queried individually against PubTator3 so that the output table contains
real literature evidence counts for every curated term rather than relying
on open-ended discovery (which missed well-known terms like SIDD/SAID/MARD
because PubTator3 only indexes titles, not abstracts).

For every seed term the script runs:
  1. Risk query     — PubTator3: "{seed_term} {phenotype}"
  2. Protective query — PubTator3: "{seed_term} {phenotype} {protective_suffix}"
     Papers unique to the protective query (absent from risk results) are
     flagged; a simple title-keyword sentiment score is computed.

Outputs (in <output_dir>/)
--------------------------
  phenotype_background.csv         — one row per seed term; all evidence columns
  phenotype_background_summary.csv — category-level aggregate counts
  phenotype_background.json        — full results dict for downstream use

The output CSV is the input to enrich_clusters_literature.py via --background_table.
Every unique term in the 'term' column will be used as a phenotype query term so
that every cluster feature (gene / SNP / HLA allele) is scored against each
background term independently.

Seeds files
-----------
  scripts/phenotype_seeds/type_2_diabetes.json  — curated T2D seeds
  scripts/phenotype_seeds/_template.json        — blank template

Usage
-----
  python scripts/phenotype_background_analysis.py \\
      --phenotype "type 2 diabetes" \\
      --output_dir results/diabetes/background \\
      [--seeds_file scripts/phenotype_seeds/type_2_diabetes.json] \\
      [--ncbi_api_key  <key>] \\
      [--open_targets]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from datetime import datetime
from typing import Any

import pandas as pd
import requests


# ── Constants ─────────────────────────────────────────────────────────────────

PUBTATOR3_BASE   = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
OPEN_TARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"

# Categories that are present in seeds files (in display order)
SEED_CATEGORIES = [
    'subtypes',
    'co_phenotypes',
    'mechanisms',
    'pathways',
    'complications',
    'comorbidities',
    'protective',
]

# Protective modifiers appended for protective queries
PROTECTIVE_SUFFIXES = [
    'protective',
    'prevention',
    'reduces risk',
    'resilience',
    'inverse association',
]

# Title-level protective sentiment keywords
PROTECTIVE_TITLE_KEYWORDS = [
    'protective', 'prevent', 'inverse', 'negat', 'reduc',
    'resilience', 'resistance', 'risk reduction', 'lower risk', 'beneficial',
]

CACHE_DIR: str | None = None


# ── HTTP helpers ───────────────────────────────────────────────────────────────

def _get(url: str, params: dict | None = None, retries: int = 3,
         wait: float = 1.5) -> Any | None:
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 200:
                try:
                    return r.json()
                except Exception:
                    return None
            if r.status_code == 429:
                time.sleep(wait * (attempt + 2))
                continue
            print(f"    [WARN] GET {url}: HTTP {r.status_code}")
            return None
        except requests.exceptions.RequestException as exc:
            print(f"    [WARN] Request attempt {attempt + 1} failed: {exc}")
            time.sleep(wait)
    return None


def _graphql_post(query: str, variables: dict) -> Any | None:
    try:
        r = requests.post(
            OPEN_TARGETS_API,
            json={'query': query, 'variables': variables},
            headers={'Content-Type': 'application/json'},
            timeout=30,
        )
        if r.status_code == 200:
            return r.json()
        print(f"    [WARN] Open Targets GraphQL: HTTP {r.status_code}")
    except requests.exceptions.RequestException as exc:
        print(f"    [WARN] Open Targets request failed: {exc}")
    return None


# ── JSON cache ─────────────────────────────────────────────────────────────────

def _cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, f'{name}.json')


def _load_cache(name: str) -> dict:
    p = _cache_path(name)
    if os.path.exists(p):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _save_cache(name: str, data: dict) -> None:
    try:
        with open(_cache_path(name), 'w') as f:
            json.dump(data, f, indent=2)
    except Exception:
        pass


# ── Seeds loader ───────────────────────────────────────────────────────────────

def _phenotype_slug(phenotype: str) -> str:
    """Convert phenotype string to filename slug (lowercase, underscores)."""
    return re.sub(r'[^a-z0-9]+', '_', phenotype.lower()).strip('_')


def load_seeds(seeds_file: str | None, phenotype: str,
               script_dir: str | None = None) -> dict:
    """
    Load a seeds JSON file.  Resolution order:
      1. Explicit --seeds_file argument
      2. Auto-detect: scripts/phenotype_seeds/<phenotype_slug>.json
      3. Empty dict with warning (graceful degradation)

    Parameters
    ----------
    seeds_file : str | None  — explicit path from CLI
    phenotype  : str         — primary phenotype string (for auto-detect)
    script_dir : str | None  — directory of this script (for auto-detect)

    Returns
    -------
    dict  keyed by category; values are lists of seed term strings
    """
    candidates = []

    if seeds_file:
        candidates.append(seeds_file)

    # Auto-detect alongside this script
    if script_dir:
        slug = _phenotype_slug(phenotype)
        candidates.append(os.path.join(script_dir, 'phenotype_seeds', f'{slug}.json'))

    for path in candidates:
        if os.path.exists(path):
            try:
                with open(path) as f:
                    raw = json.load(f)
                # Strip metadata keys; keep only category lists
                seeds = {
                    k: v for k, v in raw.items()
                    if isinstance(v, list) and k in SEED_CATEGORIES
                }
                total = sum(len(v) for v in seeds.values())
                print(f"  Seeds file        : {path}")
                print(f"  Seed terms loaded : {total} across {len(seeds)} categories")
                return seeds
            except Exception as exc:
                print(f"  [WARN] Could not parse seeds file '{path}': {exc}")

    print(f"  [WARN] No seeds file found for '{phenotype}'. "
          f"Create scripts/phenotype_seeds/{_phenotype_slug(phenotype)}.json "
          f"using scripts/phenotype_seeds/_template.json as a starting point.")
    return {}


# ── PubTator3 per-term queries ─────────────────────────────────────────────────

def _title_sentiment(title: str) -> bool:
    """Return True if the title contains a protective-sentiment keyword."""
    t = title.lower()
    return any(kw in t for kw in PROTECTIVE_TITLE_KEYWORDS)


def query_seed_term_risk(seed_term: str, phenotype: str) -> dict:
    """
    Query PubTator3 for co-occurrence of seed_term and phenotype.

    Returns
    -------
    dict  {'count', 'pmids' (set), 'sample_titles' (list[str])}
    """
    cache     = _load_cache('seed_risk')
    cache_key = f'{seed_term}__{phenotype}'
    if cache_key in cache:
        rec = cache[cache_key]
        return {**rec, 'pmids': set(rec.get('pmids', []))}

    query = f"{seed_term} {phenotype}"
    data  = _get(f"{PUBTATOR3_BASE}/search/", params={'text': query, 'page': 1})
    time.sleep(0.4)

    pmids  : set       = set()
    titles : list[str] = []
    count  = 0

    if data:
        count = data.get('count', 0)
        for r in data.get('results', []):
            pmid  = str(r.get('pmid', ''))
            title = r.get('title', '')
            if pmid:
                pmids.add(pmid)
                titles.append(title)

    rec = {'count': count, 'pmids': list(pmids), 'sample_titles': titles[:5]}
    cache[cache_key] = rec
    _save_cache('seed_risk', cache)
    return {**rec, 'pmids': pmids}


def query_seed_term_protective(seed_term: str, phenotype: str,
                                risk_pmids: set) -> dict:
    """
    Query PubTator3 for protective/prevention associations of seed_term × phenotype.

    Tries each PROTECTIVE_SUFFIX in turn and aggregates unique PMIDs.
    Papers absent from risk_pmids are flagged as unique protective findings.

    Returns
    -------
    dict  {
      'count':                int   — total unique protective PMIDs found
      'unique_count':         int   — papers not in risk query
      'sentiment_score':      float — fraction of titles with protective keywords
      'sample_prot_titles':   list  — titles with protective-sentiment keywords
      'sample_uniq_titles':   list  — titles unique to protective query
    }
    """
    cache     = _load_cache('seed_protective')
    cache_key = f'{seed_term}__{phenotype}'
    if cache_key in cache:
        return cache[cache_key]

    all_pmids  : set       = set()
    pmid_title : dict      = {}   # pmid → title
    prot_titles: list[str] = []

    for suffix in PROTECTIVE_SUFFIXES[:3]:
        query = f"{seed_term} {phenotype} {suffix}"
        data  = _get(f"{PUBTATOR3_BASE}/search/", params={'text': query, 'page': 1})
        time.sleep(0.4)
        if not data:
            continue
        for r in data.get('results', []):
            pmid  = str(r.get('pmid', ''))
            title = r.get('title', '')
            if pmid and pmid not in all_pmids:
                all_pmids.add(pmid)
                pmid_title[pmid] = title
                if _title_sentiment(title):
                    prot_titles.append(title)

    unique_pmids  = all_pmids - risk_pmids
    unique_titles = [pmid_title[p] for p in sorted(unique_pmids)
                     if p in pmid_title]
    all_titles    = list(pmid_title.values())

    rec = {
        'count':              len(all_pmids),
        'unique_count':       len(unique_pmids),
        'sentiment_score':    round(len(prot_titles) / max(1, len(all_titles)), 3),
        'sample_prot_titles': prot_titles[:5],
        'sample_uniq_titles': unique_titles[:5],
    }
    cache[cache_key] = rec
    _save_cache('seed_protective', cache)
    return rec


# ── Open Targets — phenotype-level top genes ───────────────────────────────────

def query_open_targets_top_genes(phenotype: str, max_genes: int = 30) -> list[dict]:
    """
    Fetch top genes associated with phenotype via Open Targets disease search.

    Returns list of {'gene', 'biotype', 'ot_score'} dicts.
    """
    cache = _load_cache('bg_ot_genes')
    if phenotype in cache:
        return cache[phenotype]

    search_query = """
    query SearchDisease($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }
    """
    result = _graphql_post(search_query, {'q': phenotype})
    time.sleep(0.4)
    if not result:
        return []

    hits = result.get('data', {}).get('search', {}).get('hits', [])
    if not hits:
        print(f"    [OT] No disease found for '{phenotype}'")
        return []

    pheno_lower = phenotype.lower()
    best_hit    = next(
        (h for h in hits if pheno_lower in h.get('name', '').lower()), hits[0]
    )
    disease_id   = best_hit['id']
    disease_name = best_hit.get('name', disease_id)
    print(f"    [OT] Resolved '{phenotype}' → {disease_name} ({disease_id})")

    assoc_query = """
    query DiseaseAssociations($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {index: 0, size: $size}) {
          rows {
            target { approvedSymbol biotype }
            score
          }
        }
      }
    }
    """
    assoc = _graphql_post(assoc_query, {'efoId': disease_id, 'size': max_genes})
    time.sleep(0.4)
    if not assoc:
        return []

    rows  = (assoc.get('data', {})
                  .get('disease', {})
                  .get('associatedTargets', {})
                  .get('rows', []))
    genes = [
        {
            'gene':     r['target']['approvedSymbol'],
            'biotype':  r['target'].get('biotype', ''),
            'ot_score': round(r['score'], 4),
        }
        for r in rows
    ]
    cache[phenotype] = genes
    _save_cache('bg_ot_genes', cache)
    return genes


# ── Main ───────────────────────────────────────────────────────────────────────

def run_background_analysis(
    phenotype: str,
    output_dir: str,
    seeds: dict,
    ncbi_api_key: str | None = None,
    include_open_targets: bool = True,
) -> dict:
    """
    Run seed-term background analysis for a phenotype.

    For every term in every category of the seeds dict, queries PubTator3
    for risk and protective associations individually, then writes:
      - phenotype_background.csv  (one row per seed term — all evidence columns)
      - phenotype_background_summary.csv  (one row per category)
      - phenotype_background.json

    Parameters
    ----------
    phenotype             : primary phenotype string
    output_dir            : directory for outputs
    seeds                 : dict {category: [term, ...]} from load_seeds()
    ncbi_api_key          : NCBI Entrez API key (unused here, kept for compatibility)
    include_open_targets  : query Open Targets for top disease-associated genes
    """
    global CACHE_DIR
    CACHE_DIR = os.path.join(output_dir, '.cache')
    os.makedirs(output_dir, exist_ok=True)

    total_terms = sum(len(v) for v in seeds.values())

    print("=" * 60)
    print("PHENOTYPE BACKGROUND ANALYSIS  (seed-term mode)")
    print("=" * 60)
    print(f"  Primary phenotype : {phenotype}")
    print(f"  Total seed terms  : {total_terms}")
    print(f"  Categories        : {', '.join(seeds.keys())}")
    print(f"  Output            : {output_dir}")
    print()

    entity_rows : list = []
    category_summary: dict = {}

    # ── Per-category, per-term loop ─────────────────────────────────────────
    for category in SEED_CATEGORIES:
        terms = seeds.get(category, [])
        if not terms:
            continue

        print(f"\n[{category}]  ({len(terms)} terms)")

        cat_risk_total = 0
        cat_prot_total = 0
        cat_uniq_total = 0

        for term in terms:
            risk_res = query_seed_term_risk(term, phenotype)
            prot_res = query_seed_term_protective(term, phenotype,
                                                  risk_res['pmids'])

            risk_count = risk_res['count']
            prot_count = prot_res['count']
            uniq_count = prot_res['unique_count']

            cat_risk_total += risk_count
            cat_prot_total += prot_count
            cat_uniq_total += uniq_count

            print(f"  {term:<45}  risk={risk_count:>6,}  "
                  f"prot={prot_count:>5,}  uniq={uniq_count:>4,}  "
                  f"sent={prot_res['sentiment_score']:.2f}")

            entity_rows.append({
                'category':                  category,
                'term':                      term,
                'entity_type':               'curated_seed',
                'count':                     risk_count,   # for backward compat
                'risk_count':                risk_count,
                'protective_count':          prot_count,
                'unique_protective_count':   uniq_count,
                'sentiment_score':           prot_res['sentiment_score'],
                'risk_sample_titles':        ' | '.join(
                    risk_res['sample_titles'][:3]),
                'protective_sample_titles':  ' | '.join(
                    prot_res['sample_prot_titles'][:3]),
                'unique_protective_titles':  ' | '.join(
                    prot_res['sample_uniq_titles'][:3]),
                'association_type':          'curated',
            })

        category_summary[category] = {
            'n_terms':          len(terms),
            'risk_total':       cat_risk_total,
            'protective_total': cat_prot_total,
            'unique_prot':      cat_uniq_total,
        }

    # ── Open Targets top genes ─────────────────────────────────────────────
    ot_genes: list[dict] = []
    if include_open_targets:
        print(f"\n[open_targets]")
        genes = query_open_targets_top_genes(phenotype, max_genes=30)
        ot_genes.extend(genes)
        print(f"  {phenotype}: {len(genes)} associated genes")

        for g in ot_genes[:20]:
            ot_count = round(g['ot_score'] * 100)
            entity_rows.append({
                'category':                  'open_targets_genes',
                'term':                      g['gene'],
                'entity_type':               'Gene',
                'count':                     ot_count,
                'risk_count':                ot_count,
                'protective_count':          0,
                'unique_protective_count':   0,
                'sentiment_score':           0.0,
                'risk_sample_titles':        f"OT score={g['ot_score']} "
                                             f"biotype={g['biotype']}",
                'protective_sample_titles':  '',
                'unique_protective_titles':  '',
                'association_type':          'open_targets',
            })

        category_summary['open_targets_genes'] = {
            'n_terms':          len(ot_genes),
            'risk_total':       sum(round(g['ot_score'] * 100) for g in ot_genes),
            'protective_total': 0,
            'unique_prot':      0,
        }

    # ── Save outputs ───────────────────────────────────────────────────────
    entity_df = pd.DataFrame(entity_rows)
    entity_df.sort_values(
        ['category', 'risk_count'],
        ascending=[True, False],
        inplace=True,
    )

    csv_path = os.path.join(output_dir, 'phenotype_background.csv')
    entity_df.to_csv(csv_path, index=False)
    print(f"\n  Saved phenotype_background.csv ({len(entity_df)} rows)")

    # Category summary
    summary_rows = [
        {
            'category':          cat,
            'n_terms':           v['n_terms'],
            'risk_paper_total':  v['risk_total'],
            'prot_paper_total':  v['protective_total'],
            'unique_prot_total': v['unique_prot'],
        }
        for cat, v in category_summary.items()
    ]
    summary_df  = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, 'phenotype_background_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved phenotype_background_summary.csv")

    # Full JSON
    json_path = os.path.join(output_dir, 'phenotype_background.json')
    with open(json_path, 'w') as f:
        json.dump({
            'phenotype':        phenotype,
            'run_date':         datetime.now().isoformat(),
            'category_summary': category_summary,
            'terms':            entity_rows,
            'ot_genes':         ot_genes,
        }, f, indent=2)
    print(f"  Saved phenotype_background.json")

    print(f"\nBackground analysis complete → {output_dir}/")
    return category_summary


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Pre-enrichment phenotype background analysis (seed-term mode).\n\n"
            "Queries PubTator3 for each curated seed term individually,\n"
            "producing a literature-backed evidence table with risk counts,\n"
            "protective counts, unique protective findings, and sentiment\n"
            "scores per term.\n\n"
            "Output CSV is used directly by enrich_clusters_literature.py\n"
            "via --background_table so every cluster feature is scored\n"
            "against every background term."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '--phenotype', required=True,
        help="Primary phenotype (e.g. 'type 2 diabetes')",
    )
    parser.add_argument(
        '--output_dir', required=True,
        help="Directory to write output CSVs and JSON",
    )
    parser.add_argument(
        '--seeds_file', default=None,
        help=(
            "Path to a phenotype seeds JSON file.\n"
            "If omitted, auto-detects scripts/phenotype_seeds/<phenotype_slug>.json.\n"
            "Use scripts/phenotype_seeds/_template.json as a starting point."
        ),
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help="NCBI Entrez API key (increases rate limit from 3 to 10 req/sec)",
    )
    parser.add_argument(
        '--open_targets', action='store_true',
        help=(
            "Also query Open Targets for the top phenotype-associated genes\n"
            "and add them to the background table (default: disabled)"
        ),
    )

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    seeds      = load_seeds(args.seeds_file, args.phenotype, script_dir)

    run_background_analysis(
        phenotype=args.phenotype,
        output_dir=args.output_dir,
        seeds=seeds,
        ncbi_api_key=args.ncbi_api_key,
        include_open_targets=args.open_targets,
    )
