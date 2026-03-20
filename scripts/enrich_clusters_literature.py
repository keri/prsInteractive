#!/usr/bin/env python3
"""
enrich_clusters_literature.py
==============================
Post-NMF literature enrichment for cluster features using:
  - DisGeNET REST API   — gene-disease association scores (GDA)
  - Open Targets        — integrated gene-disease evidence (GWAS, literature,
                          functional genomics); no API key required
  - PubTator3           — entity-aware publication co-occurrence (NCBI);
                          recognises genes, SNPs and diseases as biological
                          entities rather than plain keywords; no key required
  - NCBI PubMed         — keyword co-occurrence count via E-utilities
  - ClinVar             — variant clinical significance

Query strategy per gene/SNP × phenotype
----------------------------------------
1. DisGeNET GDA score (gene-disease association)
2. Open Targets association score (integrated evidence)
3. PubTator3 entity-annotated co-occurrence count
4. PubMed keyword co-occurrence count (fallback / complement)
5. ClinVar pathogenicity (SNPs and short gene symbols only)
6. Combined score = weighted sum of the above

Usage
-----
  python enrich_clusters_literature.py \\
    --cluster_dir /path/to/cohortBnmf/cardio \\
    --phenotypes "type 2 diabetes,epilepsy,cardio" \\
    --ncbi_api_key YOUR_NCBI_KEY \\
    [--disgenet_api_key YOUR_DISGENET_KEY]
"""

import argparse
import json
import math
import os
import time
from collections import defaultdict

import pandas as pd
import requests


# ── Constants ──────────────────────────────────────────────────────────────────

NCBI_BASE          = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DISGENET_BASE      = "https://www.disgenet.org/api"
OPEN_TARGETS_API   = "https://api.platform.opentargets.org/api/v4/graphql"
PUBTATOR3_BASE     = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"

# Significance threshold: if combined evidence below this, still include
# but flag as low-confidence.
MIN_SIGNIFICANT_SCORE = 0.05

# Weights for combined evidence score (must sum to 1.0)
WEIGHT_DISGENET      = 0.25
WEIGHT_OPEN_TARGETS  = 0.30
WEIGHT_PUBTATOR3     = 0.25
WEIGHT_CLINVAR       = 0.20

# ── Population context ────────────────────────────────────────────────────────
# When running enrichment on high-risk case clusters, queries are already
# risk-oriented by default.  For low-risk control clusters, protective modifiers
# are appended to PubTator3 and PubMed queries to surface resilience/prevention
# literature that would otherwise be drowned out by disease-association results.

PROTECTIVE_QUERY_SUFFIXES = [
    'protective', 'prevention', 'reduces risk', 'resilience',
    'resistance', 'inverse association', 'negatively associated',
]

# Keywords scanned in paper titles to assign a simple protective-sentiment flag.
# A paper is flagged protective if its title contains ≥1 of these phrases.
PROTECTIVE_TITLE_KEYWORDS = [
    'protective', 'prevent', 'inverse', 'negat', 'reduc', 'resilience',
    'resistance', 'risk reduction', 'lower risk', 'beneficial',
]

CACHE_DIR = None


# ── HTTP helpers ───────────────────────────────────────────────────────────────

def _get(url, params=None, headers=None, retries=3, wait=1.5):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=30)
            if r.status_code == 200:
                try:
                    return r.json()
                except Exception:
                    return r.text
            if r.status_code in (429, 503):
                time.sleep(wait * (attempt + 1))
                continue
            r.raise_for_status()
        except requests.RequestException as e:
            if attempt == retries - 1:
                print(f"    [WARN] GET {url}: {e}")
                return None
            time.sleep(wait)
    return None


def _graphql_post(query, variables, retries=3, wait=1.5):
    """POST a GraphQL query to Open Targets."""
    for attempt in range(retries):
        try:
            r = requests.post(
                OPEN_TARGETS_API,
                json={'query': query, 'variables': variables},
                headers={'Content-Type': 'application/json'},
                timeout=30,
            )
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 503):
                time.sleep(wait * (attempt + 1))
                continue
        except requests.RequestException as e:
            if attempt == retries - 1:
                print(f"    [WARN] GraphQL POST: {e}")
                return None
            time.sleep(wait)
    return None


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


# ── DisGeNET ───────────────────────────────────────────────────────────────────

def query_disgenet(gene_symbol, phenotype, api_key=None):
    """
    Query DisGeNET for gene-disease associations.

    Returns
    -------
    float  GDA score in [0, 1] (0 = not found)
    str    matched disease name or empty string
    """
    cache = _load_cache('disgenet')
    key   = f'{gene_symbol}__{phenotype}'
    if key in cache:
        return cache[key]['score'], cache[key]['notes']

    headers = {'Authorization': f'Bearer {api_key}'} if api_key else {}
    headers['Accept'] = 'application/json'

    data = _get(
        f"{DISGENET_BASE}/gda/gene/{gene_symbol}",
        params={'disease': phenotype, 'format': 'json'},
        headers=headers,
    )
    time.sleep(0.3)

    score = 0.0
    notes = ''
    if isinstance(data, list) and data:
        scores = [
            float(rec.get('score', 0))
            for rec in data
            if phenotype.lower() in rec.get('diseaseName', '').lower()
        ]
        if scores:
            score = max(scores)
            top   = sorted(data, key=lambda r: float(r.get('score', 0)), reverse=True)[0]
            notes = top.get('diseaseName', '')

    cache[key] = {'score': score, 'notes': notes}
    _save_cache('disgenet', cache)
    return score, notes


# ── Open Targets ───────────────────────────────────────────────────────────────

def query_open_targets(gene_symbol, phenotype):
    """
    Query Open Targets platform for gene-disease association evidence.

    Uses the public GraphQL API (no API key required).  Steps:
      1. Resolve gene symbol → Ensembl ID via target search
      2. Fetch top associated diseases; match against phenotype string
      3. Return the best matching association score [0, 1]

    For SNP rs-IDs the function returns 0 (Open Targets is gene-centric;
    use query_disgenet or query_pubtator3 for SNP-level evidence).

    Returns
    -------
    float  association score in [0, 1]
    str    matched disease name or empty string
    """
    if gene_symbol.startswith('rs'):
        return 0.0, ''   # SNP — not handled by Open Targets platform

    cache = _load_cache('open_targets')
    key   = f'{gene_symbol}__{phenotype}'
    if key in cache:
        return cache[key]['score'], cache[key]['disease']

    # ── Step 1: resolve gene symbol to Ensembl target ID ──────────────────
    search_q = """
    query SearchTarget($q: String!) {
        search(queryString: $q, entityNames: ["target"]) {
            hits { id }
        }
    }
    """
    resp = _graphql_post(search_q, {'q': gene_symbol})
    target_id = None
    if resp:
        hits = resp.get('data', {}).get('search', {}).get('hits', [])
        if hits:
            target_id = hits[0]['id']

    if not target_id:
        cache[key] = {'score': 0.0, 'disease': ''}
        _save_cache('open_targets', cache)
        return 0.0, ''

    # ── Step 2: fetch top disease associations ─────────────────────────────
    assoc_q = """
    query TargetDiseases($targetId: String!) {
        target(ensemblId: $targetId) {
            associatedDiseases(page: {index: 0, size: 50}) {
                rows {
                    score
                    disease { name }
                }
            }
        }
    }
    """
    resp = _graphql_post(assoc_q, {'targetId': target_id})
    time.sleep(0.4)

    score       = 0.0
    disease_hit = ''
    if resp:
        rows = (resp.get('data', {})
                    .get('target', {})
                    .get('associatedDiseases', {})
                    .get('rows', []))
        pheno_lower = phenotype.lower()
        matching = [
            r for r in rows
            if pheno_lower in r.get('disease', {}).get('name', '').lower()
        ]
        if matching:
            best        = max(matching, key=lambda r: r['score'])
            score       = best['score']
            disease_hit = best['disease']['name']
        elif rows:
            # No phenotype match — return a discounted highest score as background
            best  = max(rows, key=lambda r: r['score'])
            score = best['score'] * 0.2

    cache[key] = {'score': round(score, 6), 'disease': disease_hit}
    _save_cache('open_targets', cache)
    return round(score, 6), disease_hit


# ── PubTator3 ──────────────────────────────────────────────────────────────────

def query_pubtator3(term, phenotype):
    """
    Query PubTator3 for entity-annotated publications co-mentioning
    term (gene symbol or rs-ID) and phenotype.

    PubTator3 recognises biological entities (genes, variants, diseases,
    chemicals) in PubMed/PMC text, so the co-occurrence count is more
    specific than a plain keyword search — "BMI" is recognised as a
    clinical measure, "rs12345" as a genomic variant, etc.

    Returns
    -------
    int   number of entity-annotated publications mentioning both term and phenotype
    str   detected entity type ('Gene', 'Variant', 'Chemical', 'Disease', or '')
    """
    cache = _load_cache('pubtator3')
    key   = f'{term}__{phenotype}'
    if key in cache:
        return cache[key]['count'], cache[key]['entity_type']

    query   = f"{term} {phenotype}"
    data    = _get(
        f"{PUBTATOR3_BASE}/search/",
        # Note: 'sort' parameter removed — PubTator3 API returns 400 Bad Request
        # when sort=score is included; default relevance ordering is sufficient.
        params={'text': query, 'page': 1},
    )
    time.sleep(0.35)

    count       = 0
    entity_type = ''
    if data:
        count = data.get('count', 0)
        # Try to detect the entity type of 'term' from the first result's annotations
        for result in data.get('results', [])[:3]:
            for passage in result.get('passages', []):
                for ann in passage.get('annotations', []):
                    if term.lower() in ann.get('text', '').lower():
                        entity_type = ann.get('infons', {}).get('type', '')
                        break
                if entity_type:
                    break
            if entity_type:
                break

    cache[key] = {'count': count, 'entity_type': entity_type}
    _save_cache('pubtator3', cache)
    return count, entity_type


# ── PubMed co-occurrence (keyword fallback) ────────────────────────────────────

def query_pubmed_cooccurrence(term, phenotype, api_key=None):
    """
    Count PubMed papers mentioning (term AND phenotype) via keyword search.
    Serves as a complement to PubTator3's entity-level counts.

    Returns
    -------
    int  co-occurrence count
    """
    cache = _load_cache('pubmed_cooccurrence')
    key   = f'{term}__{phenotype}'
    if key in cache:
        return cache[key]

    query  = f'("{term}"[tiab] OR "{term}"[gene]) AND ("{phenotype}"[tiab])'
    params = {'db': 'pubmed', 'term': query, 'rettype': 'count', 'retmode': 'json'}
    if api_key:
        params['api_key'] = api_key

    data  = _get(f"{NCBI_BASE}/esearch.fcgi", params=params)
    count = 0
    if data and 'esearchresult' in data:
        count = int(data['esearchresult'].get('count', 0))
    time.sleep(0.35 if not api_key else 0.12)

    cache[key] = count
    _save_cache('pubmed_cooccurrence', cache)
    return count


# ── Protective literature query ────────────────────────────────────────────────

def query_pubtator3_protective(term, phenotype, api_key=None):
    """
    Query PubTator3 for protective/prevention associations between term and phenotype.

    Runs one query per PROTECTIVE_QUERY_SUFFIXES modifier, aggregates unique paper
    counts, and scans returned titles for protective-sentiment keywords.

    Returns
    -------
    int   total protective paper count (union across modifier queries)
    float protective sentiment score 0-1 (fraction of titles with protective keywords)
    list  sample protective titles (up to 3)
    """
    cache = _load_cache('pubtator3_protective')
    key   = f'{term}__{phenotype}'
    if key in cache:
        rec = cache[key]
        return rec['count'], rec['sentiment'], rec['titles']

    seen_pmids:  set  = set()
    prot_titles: list = []
    all_titles:  list = []

    for modifier in PROTECTIVE_QUERY_SUFFIXES:
        query = f"{term} {phenotype} {modifier}"
        data  = _get(
            f"{PUBTATOR3_BASE}/search/",
            params={'text': query, 'page': 1},
        )
        time.sleep(0.35)
        if not data:
            continue
        for result in data.get('results', []):
            pmid  = str(result.get('pmid', ''))
            title = result.get('title', '')
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                all_titles.append(title)
                title_lower = title.lower()
                if any(kw in title_lower for kw in PROTECTIVE_TITLE_KEYWORDS):
                    prot_titles.append(title)

    count     = len(seen_pmids)
    sentiment = round(len(prot_titles) / max(1, len(all_titles)), 3)
    titles    = prot_titles[:3]

    cache[key] = {'count': count, 'sentiment': sentiment, 'titles': titles}
    _save_cache('pubtator3_protective', cache)
    return count, sentiment, titles


def _detect_population_context(cluster_dir):
    """
    Auto-detect population context from the cluster_dir path.
    Returns 'high_risk' or 'low_control' or 'both'.
    """
    parts = cluster_dir.replace('\\', '/').split('/')
    if 'high_risk' in parts:
        return 'high_risk'
    if 'low_control' in parts:
        return 'low_control'
    return 'both'


# ── ClinVar ────────────────────────────────────────────────────────────────────

def query_clinvar_significance(snp_or_gene, api_key=None):
    """
    Fetch clinical significance from ClinVar for a SNP rs-ID or gene.

    Returns
    -------
    str   'pathogenic' | 'likely_pathogenic' | 'uncertain' | 'benign' | 'not_found'
    """
    cache = _load_cache('clinvar')
    if snp_or_gene in cache:
        return cache[snp_or_gene]

    params = {
        'db':      'clinvar',
        'term':    snp_or_gene,
        'retmode': 'json',
        'retmax':  5,
    }
    if api_key:
        params['api_key'] = api_key

    search = _get(f"{NCBI_BASE}/esearch.fcgi", params=params)
    time.sleep(0.35 if not api_key else 0.12)

    sig = 'not_found'
    if search and search.get('esearchresult', {}).get('idlist'):
        uid  = search['esearchresult']['idlist'][0]
        summ = _get(
            f"{NCBI_BASE}/esummary.fcgi",
            params={'db': 'clinvar', 'id': uid, 'retmode': 'json'},
        )
        time.sleep(0.35 if not api_key else 0.12)
        if summ:
            doc = summ.get('result', {}).get(uid, {})
            sig = doc.get('clinical_significance', {}).get('description', 'not_found').lower()
            sig = sig.split('/')[0].strip().replace(' ', '_')

    cache[snp_or_gene] = sig
    _save_cache('clinvar', cache)
    return sig


# ── Combined evidence scoring ─────────────────────────────────────────────────

def compute_combined_score(disgenet_score, open_targets_score,
                           pubtator3_count, clinvar_sig):
    """
    Combine evidence from DisGeNET, Open Targets, PubTator3 and ClinVar
    into a single [0, 1] score.

    PubMed keyword counts are deliberately excluded from the formula here
    (they are stored in the evidence table for reference) because PubTator3
    entity counts already capture publication evidence more precisely.

    Weights: DisGeNET=0.25, Open Targets=0.30, PubTator3=0.25, ClinVar=0.20
    """
    # PubTator3: sigmoid normalisation — ~20 entity-annotated papers → 0.63
    pubtator3_norm = 1 - math.exp(-pubtator3_count / 20.0)

    clinvar_val = {
        'pathogenic':             1.0,
        'likely_pathogenic':      0.8,
        'pathogenic/likely':      0.9,
        'uncertain_significance': 0.2,
        'benign':                 0.0,
        'likely_benign':          0.05,
        'not_found':              0.0,
    }.get(clinvar_sig, 0.1)

    combined = (
        WEIGHT_DISGENET     * disgenet_score     +
        WEIGHT_OPEN_TARGETS * open_targets_score +
        WEIGHT_PUBTATOR3    * pubtator3_norm     +
        WEIGHT_CLINVAR      * clinvar_val
    )
    return round(min(combined, 1.0), 6)


# ── Feature → term extraction ─────────────────────────────────────────────────

def extract_terms_from_feature(feature_name, snp_gene_map=None):
    """
    Given a prefixed feature column name, return queryable terms.

    Examples:
      'gen_rs12345'             → ['rs12345', 'GENE_NAME'] (if snp_gene_map provided)
      'gen_rs12345_rs67890'     → ['rs12345', 'rs67890', ...]
      'clin_bmi'                → ['bmi', 'body mass index']
      '{cohort}_prs'            → []  (PRS anchors not queried individually)
    """
    if '_prs' in feature_name and not feature_name.startswith('gen_'):
        return []

    name = feature_name
    for prefix in ('gen_', 'clin_'):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break

    # Split comma-separated interaction components (e.g. "Systolic blood pressure,DRB5_9901")
    components = [c.strip() for c in name.split(',') if c.strip()]

    terms = []
    for comp in components:
        if comp.startswith('rs') or '_rs' in comp:
            # SNP or SNP interaction pair
            parts   = comp.split('_')
            rs_ids  = [p for p in parts if p.startswith('rs')]
            terms.extend(rs_ids)
            if snp_gene_map:
                for snp in rs_ids:
                    terms.extend(snp_gene_map.get(snp, []))
        elif '_' in comp and not any(p.startswith('rs') for p in comp.split('_')):
            # Imputed HLA allele: GENE_allele → primary HLA form;
            # fallback candidates are resolved in the evidence loop.
            gene   = comp.split('_')[0]
            allele = '_'.join(comp.split('_')[1:])
            terms.append(f'HLA-{gene}*{allele}')   # e.g. HLA-DRB5*9901
        else:
            # Plain text term (clinical measure, phenotype label, etc.)
            terms.append(comp.replace('_', ' '))

    return list(dict.fromkeys(terms))


def _hla_candidates(primary_term):
    """
    Given the primary HLA query term (e.g. 'HLA-DRB5*9901'), return the ordered
    list of fallback candidates used for NCBI/literature queries:
      1. HLA-GENE*allele  (HLA-DRB5*9901)
      2. HLA-GENE         (HLA-DRB5)
      3. GENE*allele      (DRB5*9901)   ← NCBI only, not GTEx
    Returns [] if the term is not in HLA-GENE*allele format.
    """
    if not primary_term.startswith('HLA-') or '*' not in primary_term:
        return []
    rest  = primary_term[4:]                        # DRB5*9901
    gene  = rest.split('*')[0]                      # DRB5
    allele = rest.split('*')[1]                     # 9901
    return [
        primary_term,                               # HLA-DRB5*9901
        f'HLA-{gene}',                              # HLA-DRB5
        f'{gene}*{allele}',                         # DRB5*9901
    ]


def _resolve_hla_term_for_literature(primary_term, phenotype, min_count=5):
    """
    Try HLA candidate terms in order and return the first whose PubTator3
    count is >= min_count, or the last candidate if none exceed the threshold.

    Returns
    -------
    str   resolved query term
    str   resolution note for logging
    """
    candidates = _hla_candidates(primary_term)
    if not candidates:
        return primary_term, ''

    for candidate in candidates:
        count, _ = query_pubtator3(candidate, phenotype)
        if count >= min_count:
            note = (f"'{primary_term}' → '{candidate}' ({count} PubTator3 hits)"
                    if candidate != primary_term else f"({count} PubTator3 hits)")
            return candidate, note

    # All candidates below threshold — use last one (most permissive)
    fallback = candidates[-1]
    return fallback, f"'{primary_term}' low coverage → fallback '{fallback}'"


# ── Background table loader ────────────────────────────────────────────────────

def load_phenotypes_from_background_table(
    background_table_path: str,
    primary_phenotype: str | None = None,
    association_type: str | None = None,
    min_count: int = 0,
) -> list[str]:
    """
    Load all unique terms from a phenotype_background.csv produced by
    phenotype_background_analysis.py and return them as a phenotype list
    for use in enrich_with_literature().

    Every gene/SNP feature will be queried against every term in this list,
    producing a full (feature × background_term) evidence matrix so that
    each cluster feature has scores for "type 2 diabetes" AND
    "beta-cell dysfunction" AND "metabolic syndrome" etc. simultaneously.

    Parameters
    ----------
    background_table_path : str  — path to phenotype_background.csv
    primary_phenotype     : str | None — if given, it is prepended to the list
                             and any duplicate removed; ensures the main phenotype
                             is always included even if not in the table
    association_type      : str | None — filter to 'risk', 'protective', or None
                             (both); default: None = include all terms
    min_count             : int — minimum count threshold to include a term;
                             0 means include all (default)

    Returns
    -------
    list[str]  deduplicated phenotype terms in discovery order
    """
    try:
        bg = pd.read_csv(background_table_path)
    except Exception as exc:
        print(f"  [WARN] Could not read background table '{background_table_path}': {exc}")
        return [primary_phenotype] if primary_phenotype else []

    required_cols = {'term', 'count'}
    missing = required_cols - set(bg.columns)
    if missing:
        print(f"  [WARN] background table missing columns: {missing}. "
              f"Falling back to primary phenotype only.")
        return [primary_phenotype] if primary_phenotype else []

    if association_type and 'association_type' in bg.columns:
        bg = bg[bg['association_type'] == association_type]

    if min_count > 0:
        bg = bg[bg['count'] >= min_count]

    # Preserve discovery order (highest count first within each category)
    terms_ordered: list[str] = []
    seen: set[str] = set()
    for term in bg.sort_values('count', ascending=False)['term']:
        t = str(term).strip()
        if t and t.lower() not in seen:
            terms_ordered.append(t)
            seen.add(t.lower())

    # Prepend primary phenotype if supplied and not already present
    if primary_phenotype:
        pp_lower = primary_phenotype.lower()
        if pp_lower not in seen:
            terms_ordered.insert(0, primary_phenotype)

    print(f"  Background table  : {background_table_path}")
    print(f"  Phenotype terms   : {len(terms_ordered)} loaded "
          f"(every feature will be scored against each term)")

    return terms_ordered


# ── Main enrichment ───────────────────────────────────────────────────────────

def enrich_with_literature(
    cluster_dir,
    phenotypes,
    background_table=None,
    disgenet_api_key=None,
    ncbi_api_key=None,
    top_n=20,
    population_context='auto',
    append_mode=True,
):
    """
    Main literature enrichment routine.

    Queries DisGeNET, Open Targets, PubTator3, PubMed and ClinVar for each
    top feature × phenotype pair identified in the bNMF cluster output.

    Every feature (gene / SNP / HLA allele) is scored against EVERY phenotype
    term, producing a full cross-product evidence matrix.  When background_table
    is supplied, all terms from that table are used as phenotype terms so that
    each cluster feature gets a score for "type 2 diabetes" AND
    "beta-cell dysfunction" AND "metabolic syndrome" etc. simultaneously.

    Parameters
    ----------
    cluster_dir        : str — bNMF output directory (contains feature_loadings.csv)
    phenotypes         : list[str] — phenotype terms to query; used as the base list
                          and also as the primary_phenotype anchor when background_table
                          is supplied (so the primary phenotype is always included)
    background_table   : str | None — path to phenotype_background.csv produced by
                          phenotype_background_analysis.py.  When supplied, ALL unique
                          terms in the table are appended to the phenotypes list so
                          that every feature × every background term is queried.
                          No ranking or truncation is applied — all terms are used.
    disgenet_api_key   : str | None — DisGeNET API key (free at disgenet.org)
    ncbi_api_key       : str | None — NCBI API key (3→10 req/sec)
    top_n              : int — top features per cluster to query
    population_context : str — 'auto' | 'both' | 'high_risk' | 'low_control'
        Controls which additional queries are run:
          'high_risk'   — standard risk-oriented queries (default behaviour)
          'low_control' — adds protective/prevention queries per feature × phenotype;
                          output columns include protective_count, protective_sentiment,
                          and unique_protective_titles (findings absent from risk queries)
          'both'/'auto' — auto-detect from cluster_dir path; if not detectable, uses 'both'
    append_mode        : bool — if True and output CSVs already exist, new columns are
                          merged into them rather than overwriting (default: True)

    Returns
    -------
    dict  {'evidence_df', 'cluster_ranked_df', 'cluster_summary_df'}
    """
    global CACHE_DIR
    CACHE_DIR = os.path.join(cluster_dir, 'enrichment', '.cache')
    out_dir   = os.path.join(cluster_dir, 'enrichment', 'literature')
    os.makedirs(out_dir, exist_ok=True)

    # ── Resolve population context ─────────────────────────────────────────
    if population_context == 'auto':
        population_context = _detect_population_context(cluster_dir)
    run_protective = (population_context == 'low_control')
    print(f"  Population context : {population_context}"
          f"{'  [protective queries enabled]' if run_protective else ''}")

    # ── Expand phenotypes from background table ────────────────────────────
    # When a background table is provided, all unique terms from it are added
    # to the phenotypes list so every feature is scored against every term.
    # No ranking or truncation — the full cross-product is always computed.
    if background_table:
        # primary_phenotype = first element of the phenotypes list acts as anchor
        primary = phenotypes[0] if phenotypes else None
        bg_terms = load_phenotypes_from_background_table(
            background_table,
            primary_phenotype=primary,
        )
        # Merge: keep original phenotypes first, then append background terms
        # that are not already present (case-insensitive dedup)
        n_before   = len(phenotypes)
        seen_lower = {p.lower() for p in phenotypes}
        extra: list[str] = []
        for t in bg_terms:
            if t.lower() not in seen_lower:
                extra.append(t)
                seen_lower.add(t.lower())
        phenotypes = list(phenotypes) + extra
        print(f"  Total query terms : {len(phenotypes)} "
              f"({n_before} explicit + {len(extra)} from background table)")

    # ── Load feature loadings ──────────────────────────────────────────────
    loadings_path = os.path.join(cluster_dir, 'feature_loadings.csv')
    if not os.path.exists(loadings_path):
        print(f"  [ERROR] feature_loadings.csv not found in {cluster_dir}")
        return {}

    H_df = pd.read_csv(loadings_path, index_col=0)

    # Load SNP→gene map from GTEx enrichment cache if available
    snp_gene_map  = {}
    snp_gene_path = os.path.join(CACHE_DIR, 'snp_gene_map.json')
    if os.path.exists(snp_gene_path):
        snp_gene_map = _load_cache('snp_gene_map')

    # ── Collect top features per cluster ──────────────────────────────────
    cluster_top: dict = {}
    for cluster in H_df.index:
        row            = H_df.loc[cluster]
        cluster_top[cluster] = row.nlargest(top_n).index.tolist()

    all_features = list({f for feats in cluster_top.values() for f in feats})
    all_terms    = {
        feat: extract_terms_from_feature(feat, snp_gene_map)
        for feat in all_features
    }

    # ── Evidence gathering ─────────────────────────────────────────────────
    unique_pairs = sorted({
        (term, pheno)
        for terms in all_terms.values()
        for term in terms
        for pheno in phenotypes
    })

    print(f"\n  Evidence sources: DisGeNET, Open Targets, PubTator3, PubMed, ClinVar")
    print(f"  Querying {len(unique_pairs)} term × phenotype pairs...")

    evidence_rows = []
    for idx, (term, pheno) in enumerate(unique_pairs):
        if idx % 20 == 0:
            print(f"    {idx}/{len(unique_pairs)} queries completed")

        # ── HLA fallback resolution ────────────────────────────────────────
        # For HLA primary terms (HLA-GENE*allele), try progressively broader
        # candidates until >= 5 PubTator3 hits are found:
        #   1. HLA-DRB5*9901  2. HLA-DRB5  3. DRB5*9901
        query_term = term
        hla_note   = ''
        if term.startswith('HLA-') and '*' in term:
            query_term, hla_note = _resolve_hla_term_for_literature(term, pheno)
            if hla_note:
                print(f"    [HLA-lit] {hla_note}")

        # DisGeNET — gene-disease association score
        disgenet_score, disgenet_disease = query_disgenet(
            query_term, pheno, api_key=disgenet_api_key
        )

        # Open Targets — integrated multi-evidence gene-disease score
        ot_score, ot_disease = query_open_targets(query_term, pheno)

        # PubTator3 — entity-annotated co-occurrence (already fetched during HLA
        # resolution if HLA; re-use from cache otherwise)
        pt3_count, pt3_entity_type = query_pubtator3(query_term, pheno)

        # PubMed — keyword co-occurrence (stored for reference)
        pubmed_count = query_pubmed_cooccurrence(query_term, pheno, api_key=ncbi_api_key)

        # ClinVar — variant/gene pathogenicity
        is_snp_or_gene = (
            query_term.startswith('rs')
            or (query_term.replace('-', '').replace('*', '').isalpha()
                and len(query_term) <= 15)
        )
        clinvar_sig = (
            query_clinvar_significance(query_term, api_key=ncbi_api_key)
            if is_snp_or_gene else 'not_found'
        )

        combined = compute_combined_score(
            disgenet_score, ot_score, pt3_count, clinvar_sig
        )

        # ── Protective queries (low_control context only) ──────────────────
        prot_count     = 0
        prot_sentiment = 0.0
        prot_titles    = []
        unique_prot    = 0
        if run_protective:
            prot_count, prot_sentiment, prot_titles = query_pubtator3_protective(
                query_term, pheno, api_key=ncbi_api_key
            )
            # Unique protective = papers with protective signal not already in
            # the standard risk query (approximated: prot_count beyond risk count)
            unique_prot = max(0, prot_count - pt3_count)

        evidence_rows.append({
            'term':                  term,
            'query_term':            query_term,
            'phenotype':             pheno,
            'population_context':    population_context,
            'disgenet_score':        round(disgenet_score, 6),
            'disgenet_disease':      disgenet_disease,
            'open_targets_score':    ot_score,
            'open_targets_disease':  ot_disease,
            'pubtator3_count':       pt3_count,
            'pubtator3_entity_type': pt3_entity_type,
            'protective_count':      prot_count,
            'protective_sentiment':  prot_sentiment,
            'unique_protective':     unique_prot,
            'protective_titles':     ' | '.join(prot_titles) if prot_titles else '',
            'pubmed_count':          pubmed_count,
            'clinvar_sig':           clinvar_sig,
            'combined_score':        combined,
        })

    evidence_df = pd.DataFrame(evidence_rows)
    evidence_df.sort_values('combined_score', ascending=False, inplace=True)
    evidence_path = os.path.join(out_dir, 'gene_evidence_table.csv')
    if append_mode and os.path.exists(evidence_path):
        existing = pd.read_csv(evidence_path)
        new_cols = [c for c in evidence_df.columns if c not in existing.columns]
        if new_cols:
            merged = existing.merge(
                evidence_df[['term', 'phenotype'] + new_cols],
                on=['term', 'phenotype'], how='left',
            )
            merged.to_csv(evidence_path, index=False)
            print(f"    Appended {len(new_cols)} new columns to gene_evidence_table.csv")
        else:
            print(f"    gene_evidence_table.csv up to date (no new columns)")
    else:
        evidence_df.to_csv(evidence_path, index=False)
        print(f"    Saved gene_evidence_table.csv ({len(evidence_df)} rows)")

    # ── Cluster-level ranking ──────────────────────────────────────────────
    ev_lookup = (
        evidence_df.set_index(['term', 'phenotype'])['combined_score'].to_dict()
    )

    cluster_ranked_rows = []
    for cluster, top_feats in cluster_top.items():
        for rank, feat in enumerate(top_feats, start=1):
            loading = float(H_df.loc[cluster, feat])
            for term in all_terms.get(feat, []):
                for pheno in phenotypes:
                    score = ev_lookup.get((term, pheno), 0.0)
                    cluster_ranked_rows.append({
                        'cluster':            cluster,
                        'feature':            feat,
                        'term':               term,
                        'phenotype':          pheno,
                        'cluster_rank':       rank,
                        'loading':            round(loading, 6),
                        'combined_score':     score,
                        'weighted_rank_score': round(loading * score, 6),
                    })

    cluster_ranked_df = pd.DataFrame(cluster_ranked_rows)
    if not cluster_ranked_df.empty:
        cluster_ranked_df.sort_values(
            ['cluster', 'weighted_rank_score'], ascending=[True, False], inplace=True
        )
        ranked_path = os.path.join(out_dir, 'cluster_evidence_ranked.csv')
        if append_mode and os.path.exists(ranked_path):
            existing_r  = pd.read_csv(ranked_path)
            new_cols_r  = [c for c in cluster_ranked_df.columns
                           if c not in existing_r.columns]
            if new_cols_r:
                key_cols = ['cluster', 'feature', 'term', 'phenotype']
                merged_r = existing_r.merge(
                    cluster_ranked_df[key_cols + new_cols_r],
                    on=key_cols, how='left',
                )
                merged_r.to_csv(ranked_path, index=False)
                print(f"    Appended {len(new_cols_r)} new columns to cluster_evidence_ranked.csv")
            else:
                print(f"    cluster_evidence_ranked.csv up to date")
        else:
            cluster_ranked_df.to_csv(ranked_path, index=False)
            print(f"    Saved cluster_evidence_ranked.csv")

    # ── Cluster summary ────────────────────────────────────────────────────
    summary_rows = []
    if not cluster_ranked_df.empty:
        for cluster, grp in cluster_ranked_df.groupby('cluster'):
            best = grp.nlargest(1, 'weighted_rank_score').iloc[0]
            summary_rows.append({
                'cluster':             cluster,
                'top_term':            best['term'],
                'top_feature':         best['feature'],
                'dominant_phenotype':  best['phenotype'],
                'max_evidence_score':  best['combined_score'],
                'weighted_rank_score': best['weighted_rank_score'],
            })

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_path = os.path.join(out_dir, 'cluster_literature_summary.csv')
        if append_mode and os.path.exists(summary_path):
            existing_s = pd.read_csv(summary_path)
            new_cols_s = [c for c in summary_df.columns
                          if c not in existing_s.columns]
            if new_cols_s:
                merged_s = existing_s.merge(
                    summary_df[['cluster'] + new_cols_s],
                    on='cluster', how='left',
                )
                merged_s.to_csv(summary_path, index=False)
                print(f"    Appended {len(new_cols_s)} new columns to cluster_literature_summary.csv")
            else:
                print(f"    cluster_literature_summary.csv up to date")
        else:
            summary_df.to_csv(summary_path, index=False)
            print(f"    Saved cluster_literature_summary.csv")
        print("\n  Cluster dominant associations:")
        for _, row in summary_df.iterrows():
            print(f"    {row['cluster']:20s} → {row['top_term']} "
                  f"× {row['dominant_phenotype']}  "
                  f"(score={row['max_evidence_score']:.3f})")

    print(f"\n  Literature enrichment complete → {out_dir}/")
    return {
        'evidence_df':        evidence_df,
        'cluster_ranked_df':  cluster_ranked_df,
        'cluster_summary_df': summary_df,
    }


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Post-NMF literature enrichment using DisGeNET, Open Targets, "
            "PubTator3, PubMed and ClinVar."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '--cluster_dir', required=True,
        help="Path to bNMF output directory containing feature_loadings.csv",
    )
    parser.add_argument(
        '--phenotypes', required=True,
        help=(
            "Comma-separated primary phenotype terms "
            "(e.g., 'type 2 diabetes').  These are always included.\n"
            "When --background_table is also supplied, ALL terms from the\n"
            "table are appended so every feature is scored against every term."
        ),
    )
    parser.add_argument(
        '--background_table', default=None,
        help=(
            "Path to phenotype_background.csv produced by\n"
            "phenotype_background_analysis.py.  When supplied, all unique\n"
            "terms in the table (subtypes, co-phenotypes, mechanisms, pathways,\n"
            "complications, comorbidities) are added to the phenotype list so\n"
            "every gene/SNP feature is scored against every background term.\n"
            "No ranking or truncation is applied."
        ),
    )
    parser.add_argument(
        '--top_n', type=int, default=20,
        help="Top features per cluster to query (default: 20)",
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help="NCBI Entrez API key — increases rate limit from 3 to 10 req/sec",
    )
    parser.add_argument(
        '--disgenet_api_key', default=None,
        help="DisGeNET API key — free registration at disgenet.org; "
             "increases rate limits and unlocks full GDA database",
    )
    parser.add_argument(
        '--population_context', default='auto',
        choices=['auto', 'both', 'high_risk', 'low_control'],
        help=(
            "Population context for query strategy (default: auto-detect from path):\n"
            "  high_risk   — standard risk-oriented queries\n"
            "  low_control — adds protective/prevention queries; extra output columns:\n"
            "                protective_count, protective_sentiment, unique_protective,\n"
            "                protective_titles\n"
            "  auto        — detect from cluster_dir path (high_risk / low_control subdir)"
        ),
    )
    parser.add_argument(
        '--no_append', action='store_true',
        help="Overwrite existing output CSVs rather than merging new columns (default: append)",
    )

    args = parser.parse_args()

    enrich_with_literature(
        cluster_dir=args.cluster_dir,
        phenotypes=[p.strip() for p in args.phenotypes.split(',')],
        background_table=args.background_table,
        disgenet_api_key=args.disgenet_api_key,
        ncbi_api_key=args.ncbi_api_key,
        top_n=args.top_n,
        population_context=args.population_context,
        append_mode=not args.no_append,
    )
