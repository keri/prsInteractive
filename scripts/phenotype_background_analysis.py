#!/usr/bin/env python3
"""
phenotype_background_analysis.py
=================================
Pre-enrichment background knowledge table for a phenotype of interest.

Before running cluster-level enrichment, this script builds a structured
literature summary covering six clinical/molecular categories for the
phenotype.  Each category is queried with:

  1. Risk query     — standard association search (phenotype + category modifier)
  2. Protective query — same search augmented with protective/prevention modifiers;
                        papers unique to this query (not found in risk results)
                        are flagged as potentially novel protective findings.

A simple title-level sentiment score (fraction of titles containing protective
keywords) distinguishes true protective associations from incidental mentions.

Categories
----------
  subtypes       — disease subgroups, endotypes, heterogeneity
  co_phenotypes  — co-occurring phenotypes, phenotypic overlap
  mechanisms     — molecular mechanisms, pathogenesis, pathophysiology
  pathways       — biological pathways, signaling, metabolic processes
  complications  — outcomes, sequelae, disease progression
  comorbidities  — comorbidities, multimorbidity

Outputs (in <output_dir>/)
--------------------------
  phenotype_background.csv   — entity-level table; one row per
                                (category, term, risk/protective)
  phenotype_background_summary.csv — category-level aggregate counts
  phenotype_background.json  — full result dict (cached for downstream use)

Usage
-----
  python scripts/phenotype_background_analysis.py \\
      --phenotype "type 2 diabetes" \\
      --output_dir results/diabetes/background \\
      [--additional_phenotypes "insulin resistance,metabolic syndrome"] \\
      [--ncbi_api_key  <key>] \\
      [--open_targets]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime
from typing import Any

import pandas as pd
import requests


# ── Constants ─────────────────────────────────────────────────────────────────

PUBTATOR3_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
NCBI_BASE      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OPEN_TARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"

# Category modifiers: each is appended to "{phenotype} {modifier}" for querying
BACKGROUND_CATEGORIES: dict[str, list[str]] = {
    'subtypes':      ['subtype', 'subgroup', 'heterogeneity', 'endotype',
                      'classification', 'stratification'],
    'co_phenotypes': ['comorbid phenotype', 'co-occurring', 'associated condition',
                      'phenotypic overlap', 'pleiotropic'],
    'mechanisms':    ['mechanism', 'pathogenesis', 'pathophysiology',
                      'molecular basis', 'etiology'],
    'pathways':      ['pathway', 'signaling pathway', 'metabolic pathway',
                      'biological process', 'gene network'],
    'complications': ['complication', 'outcome', 'sequelae', 'progression',
                      'end-organ damage'],
    'comorbidities': ['comorbidity', 'multimorbidity', 'co-morbid',
                      'concurrent disease'],
}

# Protective modifiers appended for protective queries
PROTECTIVE_SUFFIXES = [
    'protective', 'prevention', 'reduces risk',
    'resilience', 'resistance', 'inverse association',
]

# Title-level protective sentiment keywords
PROTECTIVE_TITLE_KEYWORDS = [
    'protective', 'prevent', 'inverse', 'negat', 'reduc',
    'resilience', 'resistance', 'risk reduction', 'lower risk', 'beneficial',
]

# Minimum PubTator3 count to include an entity in the output table
MIN_COUNT_THRESHOLD = 3

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
            print(f"    [WARN] Request attempt {attempt+1} failed: {exc}")
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


# ── PubTator3 queries ──────────────────────────────────────────────────────────

def _pubtator3_search(query: str, page: int = 1) -> dict | None:
    """Single PubTator3 search; returns raw API response or None."""
    return _get(
        f"{PUBTATOR3_BASE}/search/",
        params={'text': query, 'page': page},
    )


def _title_sentiment(title: str) -> bool:
    """Return True if the title contains protective-sentiment keywords."""
    t = title.lower()
    return any(kw in t for kw in PROTECTIVE_TITLE_KEYWORDS)


def _extract_entities(results: list[dict], max_entities: int = 30) -> list[dict]:
    """
    Extract annotated biological entities from PubTator3 search results.

    Returns list of {'text', 'type', 'pmid', 'title'} dicts.
    """
    entities = []
    for result in results:
        pmid  = str(result.get('pmid', ''))
        title = result.get('title', '')
        for passage in result.get('passages', []):
            for ann in passage.get('annotations', []):
                etype = ann.get('infons', {}).get('type', '')
                etext = ann.get('text', '').strip()
                if etext and etype and etype not in ('Species',):
                    entities.append({
                        'text':  etext,
                        'type':  etype,
                        'pmid':  pmid,
                        'title': title,
                    })
                if len(entities) >= max_entities:
                    return entities
    return entities


def query_category_risk(phenotype: str, category: str,
                         modifiers: list[str],
                         api_key: str | None = None) -> dict:
    """
    Query PubTator3 for risk associations of phenotype within a category.

    Returns
    -------
    dict {
      'count':     int   — total paper count
      'entities':  list  — entity dicts {text, type, pmid, title}
      'pmids':     set   — set of PMIDs found
    }
    """
    cache = _load_cache('bg_risk')
    cache_key = f'{phenotype}__{category}'
    if cache_key in cache:
        rec = cache[cache_key]
        return {**rec, 'pmids': set(rec.get('pmids', []))}

    all_entities: list[dict] = []
    all_pmids:    set  = set()
    total_count = 0

    for modifier in modifiers:
        query = f"{phenotype} {modifier}"
        data  = _pubtator3_search(query)
        time.sleep(0.4)
        if not data:
            continue
        total_count += data.get('count', 0)
        results = data.get('results', [])
        all_entities.extend(_extract_entities(results))
        for r in results:
            pmid = str(r.get('pmid', ''))
            if pmid:
                all_pmids.add(pmid)

    rec = {
        'count':    total_count,
        'entities': all_entities[:60],
        'pmids':    list(all_pmids),
    }
    cache[cache_key] = rec
    _save_cache('bg_risk', cache)
    return {**rec, 'pmids': all_pmids}


def query_category_protective(phenotype: str, category: str,
                               modifiers: list[str],
                               risk_pmids: set,
                               api_key: str | None = None) -> dict:
    """
    Query PubTator3 for protective associations within a category.

    Subtracts risk_pmids to find unique protective findings.

    Returns
    -------
    dict {
      'count':               int   — total protective paper count
      'unique_count':        int   — papers not found in risk query
      'sentiment_score':     float — fraction of titles with protective keywords
      'sample_prot_titles':  list  — up to 5 unique protective titles
      'entities':            list  — entity dicts
    }
    """
    cache = _load_cache('bg_protective')
    cache_key = f'{phenotype}__{category}'
    if cache_key in cache:
        return cache[cache_key]

    all_entities:  list[dict] = []
    all_pmids:     set  = set()
    prot_titles:   list = []
    all_titles:    list = []

    for modifier in modifiers:
        for prot_suffix in PROTECTIVE_SUFFIXES[:3]:  # top 3 protective suffixes
            query = f"{phenotype} {modifier} {prot_suffix}"
            data  = _pubtator3_search(query)
            time.sleep(0.4)
            if not data:
                continue
            for result in data.get('results', []):
                pmid  = str(result.get('pmid', ''))
                title = result.get('title', '')
                if pmid and pmid not in all_pmids:
                    all_pmids.add(pmid)
                    all_titles.append(title)
                    if _title_sentiment(title):
                        prot_titles.append(title)
                    all_entities.extend(
                        _extract_entities([result], max_entities=5)
                    )

    unique_pmids  = all_pmids - risk_pmids
    unique_titles = [t for pmid, t in zip(all_pmids, all_titles)
                     if pmid in unique_pmids]

    rec = {
        'count':              len(all_pmids),
        'unique_count':       len(unique_pmids),
        'sentiment_score':    round(len(prot_titles) / max(1, len(all_titles)), 3),
        'sample_prot_titles': prot_titles[:5],
        'sample_uniq_titles': unique_titles[:5],
        'entities':           all_entities[:60],
    }
    cache[cache_key] = rec
    _save_cache('bg_protective', cache)
    return rec


# ── Open Targets — phenotype-level top genes ───────────────────────────────────

def query_open_targets_top_genes(phenotype: str, max_genes: int = 30) -> list[dict]:
    """
    Fetch top genes associated with phenotype via Open Targets disease search.

    Returns list of {'gene', 'score', 'category'} dicts.
    """
    cache = _load_cache('bg_ot_genes')
    if phenotype in cache:
        return cache[phenotype]

    # Step 1: find disease EFO ID
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

    # Pick best hit — prefer exact-ish match
    pheno_lower = phenotype.lower()
    best_hit = next(
        (h for h in hits if pheno_lower in h.get('name', '').lower()),
        hits[0]
    )
    disease_id   = best_hit['id']
    disease_name = best_hit.get('name', disease_id)
    print(f"    [OT] Resolved '{phenotype}' → {disease_name} ({disease_id})")

    # Step 2: fetch top associated genes
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

    rows = (assoc.get('data', {})
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


# ── Entity aggregation ─────────────────────────────────────────────────────────

def _aggregate_entities(entities: list[dict]) -> list[dict]:
    """
    Deduplicate and count entities by (text, type).

    Returns list sorted by count descending.
    """
    counts: dict[tuple, int] = {}
    samples: dict[tuple, str] = {}
    for e in entities:
        key = (e['text'], e['type'])
        counts[key] = counts.get(key, 0) + 1
        samples.setdefault(key, e.get('title', ''))
    return [
        {'term': k[0], 'entity_type': k[1], 'count': v, 'sample_title': samples[k]}
        for k, v in sorted(counts.items(), key=lambda x: -x[1])
    ]


# ── Main ───────────────────────────────────────────────────────────────────────

def run_background_analysis(
    phenotype: str,
    output_dir: str,
    additional_phenotypes: list[str] | None = None,
    ncbi_api_key: str | None = None,
    include_open_targets: bool = True,
) -> dict:
    """
    Run the full pre-enrichment background analysis for a phenotype.

    Parameters
    ----------
    phenotype             : primary phenotype string (e.g. 'type 2 diabetes')
    output_dir            : directory for output CSVs and JSON
    additional_phenotypes : optional list of related phenotype aliases to include
    ncbi_api_key          : NCBI Entrez API key (increases rate limit)
    include_open_targets  : if True, also query Open Targets for top associated genes

    Returns
    -------
    dict  full results keyed by category
    """
    global CACHE_DIR
    CACHE_DIR = os.path.join(output_dir, '.cache')
    os.makedirs(output_dir, exist_ok=True)

    all_phenotypes = [phenotype] + (additional_phenotypes or [])

    print("=" * 60)
    print("PHENOTYPE BACKGROUND ANALYSIS")
    print("=" * 60)
    print(f"  Primary phenotype    : {phenotype}")
    if additional_phenotypes:
        print(f"  Additional aliases   : {', '.join(additional_phenotypes)}")
    print(f"  Categories           : {', '.join(BACKGROUND_CATEGORIES)}")
    print(f"  Output               : {output_dir}")
    print()

    results:     dict = {}
    entity_rows: list = []

    # ── Category loop ──────────────────────────────────────────────────────
    for category, modifiers in BACKGROUND_CATEGORIES.items():
        print(f"\n[{category}]")

        cat_risk_entities:  list = []
        cat_prot_entities:  list = []
        cat_risk_count     = 0
        cat_prot_count     = 0
        cat_unique_prot    = 0
        cat_risk_pmids:    set  = set()
        cat_prot_sentiment = 0.0
        cat_prot_titles:   list = []
        cat_uniq_titles:   list = []

        for pheno in all_phenotypes:
            # Risk
            risk_res = query_category_risk(pheno, category, modifiers, ncbi_api_key)
            cat_risk_count  += risk_res['count']
            cat_risk_pmids  |= risk_res['pmids']
            cat_risk_entities.extend(risk_res['entities'])
            print(f"  risk   [{pheno}]: {risk_res['count']:,} papers  "
                  f"({len(risk_res['pmids'])} unique PMIDs)")

            # Protective
            prot_res = query_category_protective(
                pheno, category, modifiers, risk_res['pmids'], ncbi_api_key
            )
            cat_prot_count   += prot_res['count']
            cat_unique_prot  += prot_res['unique_count']
            cat_prot_sentiment = round(
                (cat_prot_sentiment + prot_res['sentiment_score']) / 2, 3
            )
            cat_prot_entities.extend(prot_res['entities'])
            cat_prot_titles.extend(prot_res['sample_prot_titles'])
            cat_uniq_titles.extend(prot_res['sample_uniq_titles'])
            print(f"  prot   [{pheno}]: {prot_res['count']:,} papers  "
                  f"unique={prot_res['unique_count']}  "
                  f"sentiment={prot_res['sentiment_score']:.2f}")

        # Aggregate entities for this category
        agg_risk = _aggregate_entities(cat_risk_entities)
        agg_prot = _aggregate_entities(cat_prot_entities)

        # Build entity rows for output table
        seen_terms: set = set()
        for ent in agg_risk:
            if ent['count'] < MIN_COUNT_THRESHOLD:
                continue
            term = ent['term']
            seen_terms.add(term)
            entity_rows.append({
                'category':            category,
                'term':                term,
                'entity_type':         ent['entity_type'],
                'association_type':    'risk',
                'count':               ent['count'],
                'sample_title':        ent['sample_title'][:200],
            })

        for ent in agg_prot:
            if ent['count'] < MIN_COUNT_THRESHOLD:
                continue
            entity_rows.append({
                'category':            category,
                'term':                ent['term'],
                'entity_type':         ent['entity_type'],
                'association_type':    'protective',
                'count':               ent['count'],
                'sample_title':        ent['sample_title'][:200],
            })

        results[category] = {
            'risk_count':           cat_risk_count,
            'protective_count':     cat_prot_count,
            'unique_protective':    cat_unique_prot,
            'sentiment_score':      cat_prot_sentiment,
            'sample_prot_titles':   cat_prot_titles[:5],
            'sample_uniq_titles':   cat_uniq_titles[:5],
            'top_risk_entities':    agg_risk[:20],
            'top_prot_entities':    agg_prot[:20],
        }

    # ── Open Targets top genes ─────────────────────────────────────────────
    ot_genes: list[dict] = []
    if include_open_targets:
        print(f"\n[open_targets]")
        for pheno in all_phenotypes:
            genes = query_open_targets_top_genes(pheno, max_genes=30)
            ot_genes.extend(genes)
            print(f"  {pheno}: {len(genes)} associated genes")

        # Deduplicate by gene symbol, keep highest score
        seen_genes: dict[str, dict] = {}
        for g in ot_genes:
            sym = g['gene']
            if sym not in seen_genes or g['ot_score'] > seen_genes[sym]['ot_score']:
                seen_genes[sym] = g
        ot_genes = sorted(seen_genes.values(), key=lambda x: -x['ot_score'])

        # Add OT genes as entity rows in 'mechanisms' category
        for g in ot_genes[:20]:
            entity_rows.append({
                'category':         'open_targets_genes',
                'term':             g['gene'],
                'entity_type':      'Gene',
                'association_type': 'risk',
                'count':            round(g['ot_score'] * 100),
                'sample_title':     f"OT score={g['ot_score']} biotype={g['biotype']}",
            })

        results['open_targets_genes'] = {
            'genes': ot_genes,
            'n':     len(ot_genes),
        }

    # ── Save outputs ───────────────────────────────────────────────────────
    entity_df = pd.DataFrame(entity_rows)
    entity_df.sort_values(['category', 'association_type', 'count'],
                          ascending=[True, True, False], inplace=True)
    entity_path = os.path.join(output_dir, 'phenotype_background.csv')
    entity_df.to_csv(entity_path, index=False)
    print(f"\n  Saved phenotype_background.csv ({len(entity_df)} rows)")

    # Category-level summary
    summary_rows = []
    for cat, r in results.items():
        if cat == 'open_targets_genes':
            continue
        summary_rows.append({
            'category':          cat,
            'risk_paper_count':  r['risk_count'],
            'prot_paper_count':  r['protective_count'],
            'unique_prot_count': r['unique_protective'],
            'prot_sentiment':    r['sentiment_score'],
            'top_risk_terms':    ', '.join(
                e['term'] for e in r['top_risk_entities'][:5]
            ),
            'top_prot_terms':    ', '.join(
                e['term'] for e in r['top_prot_entities'][:5]
            ),
            'sample_prot_title': r['sample_prot_titles'][0]
                                 if r['sample_prot_titles'] else '',
            'sample_uniq_prot':  r['sample_uniq_titles'][0]
                                 if r['sample_uniq_titles'] else '',
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, 'phenotype_background_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved phenotype_background_summary.csv")

    # Full JSON
    json_path = os.path.join(output_dir, 'phenotype_background.json')
    with open(json_path, 'w') as f:
        # Convert sets to lists for JSON serialisation
        json_safe = {}
        for cat, r in results.items():
            json_safe[cat] = {
                k: list(v) if isinstance(v, set) else v
                for k, v in r.items()
            }
        json.dump({
            'phenotype':   phenotype,
            'run_date':    datetime.now().isoformat(),
            'categories':  json_safe,
        }, f, indent=2)
    print(f"  Saved phenotype_background.json")

    print(f"\nBackground analysis complete → {output_dir}/")
    return results


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Pre-enrichment phenotype background analysis.\n\n"
            "Builds a literature-backed table of subtypes, co-phenotypes,\n"
            "molecular mechanisms, pathways, complications and comorbidities\n"
            "for a phenotype of interest.  Risk AND protective associations\n"
            "are queried separately; protective-unique findings are flagged."
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
        '--additional_phenotypes', default=None,
        help="Comma-separated phenotype aliases to include alongside the primary\n"
             "(e.g. 'T2D,diabetes mellitus type 2,insulin resistance')",
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help="NCBI Entrez API key — increases rate limit from 3 to 10 req/sec",
    )
    parser.add_argument(
        '--open_targets', action='store_true',
        help="Also query Open Targets for top phenotype-associated genes\n"
             "(adds 'open_targets_genes' section to output; default: disabled)",
    )

    args = parser.parse_args()

    extra = (
        [p.strip() for p in args.additional_phenotypes.split(',')]
        if args.additional_phenotypes else None
    )

    run_background_analysis(
        phenotype=args.phenotype,
        output_dir=args.output_dir,
        additional_phenotypes=extra,
        ncbi_api_key=args.ncbi_api_key,
        include_open_targets=args.open_targets,
    )
