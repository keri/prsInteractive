#!/usr/bin/env python3
"""
build_seeds_from_articles.py
============================
Build a phenotype seeds JSON file from research articles in a Zotero collection
or a plain-text PMID list.

Steps
-----
1. Load paper metadata (title, abstract, PMID) from Zotero SQLite or PMID file.
2. Fetch entity annotations from PubTator3 biocjson endpoint for papers with PMIDs.
3. Use Claude API (claude-opus-4-6) to classify extracted entities into the seven
   seed categories (subtypes, co_phenotypes, mechanisms, pathways, complications,
   comorbidities, protective).
4. Write a seeds JSON file consumable by phenotype_background_analysis.py.

Collection-to-category mapping
-------------------------------
If the target Zotero collection has sub-collections whose names match category
keywords (e.g. "Subtypes", "Mechanisms", "Protective factors"), papers in each
sub-collection are mapped to that category before classification — giving Claude
a strong prior that improves accuracy.

If the collection is flat (no sub-collections), all papers are pooled and Claude
classifies each entity across all categories.

RAM note: With 128 GB RAM, thousands of papers can be held in memory simultaneously.
The practical bottleneck is PubTator3 API rate limits (~3 req/s without NCBI key,
~10 req/s with key) and Claude API throughput.

Usage
-----
  # From a Zotero collection (auto-detects ~/Zotero/zotero.sqlite)
  python scripts/build_seeds_from_articles.py \\
      --phenotype "type 2 diabetes" \\
      --zotero_collection "T2D Research" \\
      --output scripts/phenotype_seeds/type_2_diabetes.json

  # Map sub-collections to specific categories
  python scripts/build_seeds_from_articles.py \\
      --phenotype "type 2 diabetes" \\
      --zotero_collection "T2D Research" \\
      --category_map "T2D Subtypes:subtypes,T2D Mechanisms:mechanisms" \\
      --output scripts/phenotype_seeds/type_2_diabetes.json

  # From a plain PMID file (one PMID per line)
  python scripts/build_seeds_from_articles.py \\
      --phenotype "type 2 diabetes" \\
      --pmids_file pmids.txt \\
      --output scripts/phenotype_seeds/type_2_diabetes.json

  # List all Zotero collections (no other args needed)
  python scripts/build_seeds_from_articles.py --list_collections

  # Merge into existing seeds file rather than overwrite
  python scripts/build_seeds_from_articles.py \\
      --phenotype "type 2 diabetes" \\
      --zotero_collection "T2D Research" \\
      --output scripts/phenotype_seeds/type_2_diabetes.json \\
      --merge
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import requests

# ── Optional Claude API ────────────────────────────────────────────────────────
try:
    import anthropic as _anthropic
    _ANTHROPIC_AVAILABLE = True
except ImportError:
    _ANTHROPIC_AVAILABLE = False

# ── Constants ─────────────────────────────────────────────────────────────────

PUBTATOR3_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"

DEFAULT_ZOTERO_DB = os.path.expanduser("~/Zotero/zotero.sqlite")

# Default seed categories — used as a fallback when no Zotero sub-folders are
# found (flat collection or PMID file).  When sub-folders ARE present, the
# actual category list is derived dynamically from the folder names so any
# phenotype with any folder structure works without changing this constant.
SEED_CATEGORIES = [
    'subtypes',       # disease subtypes and endotypes (SIDD, SIRD, MARD …)
    'mechanisms',     # molecular / cellular pathogenesis
    'pathways',       # biological signalling and metabolic pathways
    'pharma',         # pharmaceutical / drug-based interventions
    'protective',     # lifestyle protective factors (diet, exercise …)
    'comorbidities',  # co-occurring diseases, complications, comorbid conditions
]

# Human-readable descriptions used in the Claude classification prompt.
# Sub-folder name → category slug is DIRECT (no keyword heuristics);
# these descriptions are only for guiding Claude's entity assignment.
CATEGORY_DESCRIPTIONS: dict[str, str] = {
    'subtypes': (
        "disease subtypes, endotypes, subgroups, and classification systems "
        "(e.g. SIDD, SIRD, MARD, MOD, MODY, LADA, ketosis-prone diabetes)"
    ),
    'mechanisms': (
        "molecular/cellular mechanisms and pathogenesis "
        "(e.g. insulin resistance, beta-cell dysfunction, glucotoxicity, "
        "oxidative stress, ER stress, mitochondrial dysfunction)"
    ),
    'pathways': (
        "biological signalling and metabolic pathways "
        "(e.g. PI3K/AKT, AMPK, mTOR, NF-κB, gluconeogenesis, FOXO1, "
        "ceramide pathway, unfolded protein response)"
    ),
    'pharma': (
        "pharmaceutical and drug-based interventions — prescription medications "
        "and surgical procedures used in disease management "
        "(e.g. metformin, GLP-1 receptor agonists, SGLT2 inhibitors, "
        "DPP-4 inhibitors, insulin, thiazolidinediones, bariatric surgery)"
    ),
    'protective': (
        "lifestyle and non-pharmacological protective factors "
        "(e.g. Mediterranean diet, physical activity, caloric restriction, "
        "intermittent fasting, plant-based diet, resistance training, "
        "dietary fibre, omega-3 fatty acids, vitamin D) — "
        "EXCLUDE drugs (→ pharma) and clinical measures (→ comorbidities)"
    ),
    'comorbidities': (
        "co-occurring diseases, disease complications, comorbid conditions, "
        "and clinically assessed outcomes "
        "(e.g. obesity, hypertension, dyslipidaemia, NAFLD, diabetic "
        "nephropathy, retinopathy, neuropathy, cardiovascular disease, "
        "heart failure, stroke, sleep disorder, depression, CKD, cancer)"
    ),
}

# Entity types to exclude from seeds (too generic or non-biological)
EXCLUDED_ENTITY_TYPES = {'Species', 'CellLine', 'Strain'}

# Minimum entity frequency (across papers) to include in classification prompt
MIN_ENTITY_FREQ = 2

# Batch size for PubTator3 biocjson requests
PUBTATOR3_BATCH_SIZE = 50


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


# ── Zotero SQLite reader ───────────────────────────────────────────────────────

def _open_zotero_db(db_path: str) -> sqlite3.Connection:
    """
    Open Zotero SQLite in read-only mode (works even when Zotero is running).
    Uses immutable=1 to bypass file locks held by the Zotero application.
    """
    uri = f"file:{db_path}?mode=ro&immutable=1"
    try:
        conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.OperationalError as exc:
        raise RuntimeError(
            f"Cannot open Zotero database at '{db_path}'.\n"
            f"If Zotero is open, try closing it first, or this is a permissions issue.\n"
            f"Error: {exc}"
        )


def list_collections(db_path: str) -> list[dict]:
    """
    Return all Zotero collections with their full path and item counts.

    Each dict contains: collectionID, collectionName, parentCollectionID,
    full_path (e.g. 'seedCollection/T2D/Subtypes'), item_count.
    """
    conn = _open_zotero_db(db_path)
    try:
        rows = conn.execute("""
            SELECT
                c.collectionID,
                c.collectionName,
                c.parentCollectionID,
                COUNT(ci.itemID) AS item_count
            FROM collections c
            LEFT JOIN collectionItems ci ON c.collectionID = ci.collectionID
            GROUP BY c.collectionID
            ORDER BY c.collectionName
        """).fetchall()
        records = [dict(r) for r in rows]

        # Build full paths using a lookup map
        id_to_name   = {r['collectionID']: r['collectionName'] for r in records}
        id_to_parent = {r['collectionID']: r['parentCollectionID'] for r in records}

        def _full_path(cid: int) -> str:
            parts: list[str] = []
            while cid is not None:
                parts.append(id_to_name.get(cid, '?'))
                cid = id_to_parent.get(cid)
            return '/'.join(reversed(parts))

        for r in records:
            r['full_path'] = _full_path(r['collectionID'])

        # Sort by full path for a clean tree-like display
        records.sort(key=lambda r: r['full_path'].lower())
        return records
    finally:
        conn.close()


def _resolve_collection_path(conn: sqlite3.Connection,
                              path: str) -> list[int]:
    """
    Resolve a collection path to a list of matching collectionIDs.

    Supports two forms:
      - Simple name : 'T2D'
          → all collections named T2D (any level)
      - Slash path  : 'seedCollection/T2D'
          → only the T2D whose immediate parent is seedCollection

    Paths may have arbitrary depth: 'grandparent/parent/child'.
    Returns [] if no match is found.
    """
    parts = [p.strip() for p in path.split('/') if p.strip()]
    if not parts:
        return []

    if len(parts) == 1:
        # Simple name: exact match first, partial fallback
        rows = conn.execute(
            "SELECT collectionID FROM collections "
            "WHERE LOWER(collectionName) = LOWER(?)",
            (parts[0],),
        ).fetchall()
        if not rows:
            rows = conn.execute(
                "SELECT collectionID FROM collections "
                "WHERE LOWER(collectionName) LIKE LOWER(?)",
                (f'%{parts[0]}%',),
            ).fetchall()
        return [r[0] for r in rows]

    # Walk the path segment by segment
    # Start: root-level collections (no parent) matching parts[0]
    current_ids: list[int] = [
        r[0] for r in conn.execute(
            "SELECT collectionID FROM collections "
            "WHERE LOWER(collectionName) = LOWER(?) "
            "AND parentCollectionID IS NULL",
            (parts[0],),
        ).fetchall()
    ]
    if not current_ids:
        return []

    for part in parts[1:]:
        if not current_ids:
            return []
        placeholders = ','.join('?' * len(current_ids))
        current_ids = [
            r[0] for r in conn.execute(
                f"SELECT collectionID FROM collections "
                f"WHERE LOWER(collectionName) = LOWER(?) "
                f"AND parentCollectionID IN ({placeholders})",
                [part] + current_ids,
            ).fetchall()
        ]
    return current_ids


def _expand_to_descendants(conn: sqlite3.Connection,
                            root_ids: list[int]) -> list[int]:
    """Return root_ids + all descendant collectionIDs via recursive CTE."""
    if not root_ids:
        return []
    placeholders = ','.join('?' * len(root_ids))
    rows = conn.execute(f"""
        WITH RECURSIVE tree(collectionID) AS (
            SELECT collectionID FROM collections
            WHERE collectionID IN ({placeholders})
            UNION ALL
            SELECT c.collectionID FROM collections c
            JOIN tree t ON c.parentCollectionID = t.collectionID
        )
        SELECT collectionID FROM tree
    """, root_ids).fetchall()
    return [r[0] for r in rows]


def _get_subcollections_by_ids(conn: sqlite3.Connection,
                                parent_ids: list[int]) -> list[dict]:
    """Return immediate sub-collections for a list of parent collectionIDs."""
    if not parent_ids:
        return []
    placeholders = ','.join('?' * len(parent_ids))
    rows = conn.execute(
        f"SELECT collectionID, collectionName "
        f"FROM collections "
        f"WHERE parentCollectionID IN ({placeholders})",
        parent_ids,
    ).fetchall()
    return [dict(r) for r in rows]


def _query_papers_in_collections(conn: sqlite3.Connection,
                                  collection_ids: list[int]) -> list[dict]:
    """
    Fetch paper metadata (title, abstract, PMID, DOI) for items in the given
    collection IDs.  Notes and attachments are excluded.
    """
    if not collection_ids:
        return []

    placeholders = ','.join('?' * len(collection_ids))
    rows = conn.execute(f"""
        SELECT DISTINCT
            i.itemID,
            it.typeName                                                AS item_type,
            MAX(CASE WHEN f.fieldName = 'title'        THEN idv.value END) AS title,
            MAX(CASE WHEN f.fieldName = 'abstractNote' THEN idv.value END) AS abstract,
            MAX(CASE WHEN f.fieldName = 'DOI'          THEN idv.value END) AS doi,
            MAX(CASE WHEN f.fieldName = 'extra'        THEN idv.value END) AS extra,
            MAX(CASE WHEN f.fieldName = 'url'          THEN idv.value END) AS url
        FROM items i
        JOIN itemTypes it  ON i.itemTypeID = it.itemTypeID
        JOIN collectionItems ci ON i.itemID = ci.itemID
        LEFT JOIN itemData id    ON i.itemID  = id.itemID
        LEFT JOIN fields f       ON id.fieldID = f.fieldID
        LEFT JOIN itemDataValues idv ON id.valueID = idv.valueID
        WHERE ci.collectionID IN ({placeholders})
          AND it.typeName NOT IN ('note', 'attachment')
        GROUP BY i.itemID
    """, collection_ids).fetchall()

    papers = []
    for r in rows:
        pmid = _extract_pmid(r['extra'] or '', r['url'] or '')
        papers.append({
            'itemID':   r['itemID'],
            'title':    r['title'] or '',
            'abstract': r['abstract'] or '',
            'doi':      r['doi'] or '',
            'pmid':     pmid,
        })
    return papers


def _extract_pmid(extra: str, url: str) -> str:
    """
    Parse PMID from Zotero 'extra' field (e.g. 'PMID: 12345678')
    or from a PubMed URL.
    """
    if extra:
        m = re.search(r'PMID:\s*(\d{6,})', extra, re.IGNORECASE)
        if m:
            return m.group(1)
    if url:
        m = re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,})', url)
        if m:
            return m.group(1)
    return ''


def _folder_name_to_category_slug(name: str) -> str:
    """
    Convert a Zotero sub-folder name to a seed category slug.

    Only typo corrections and genuine abbreviations are normalised here.
    Semantic remapping (e.g. 'Clinical' → 'comorbidities') is intentionally
    NOT done — each folder becomes its own category so the script scales to
    any phenotype with any folder structure.

    Examples:
      'Subtypes'       → 'subtypes'
      'Pharma'         → 'pharma'
      'Pharmaceutical' → 'pharma'      (abbreviation)
      'Protective'     → 'protective'
      'Lifestyle'      → 'protective'  (synonym)
      'Mechanisms'     → 'mechanisms'
      'Mechanism'      → 'mechanisms'  (singular → plural)
      'Pathways'       → 'pathways'
      'Comorbidities'  → 'comorbidities'
      'Comorbities'    → 'comorbidities'  (typo)
      'Clinical'       → 'clinical'    (own category — NOT remapped)
      'Complications'  → 'complications'  (own category — NOT remapped)
    """
    slug = re.sub(r'[^a-z0-9]+', '_', name.strip().lower()).strip('_')
    # Only typos and clear abbreviations — no semantic remapping
    _ALIASES: dict[str, str] = {
        # pharma abbreviations
        'pharmaceutical':  'pharma',
        'pharmaceuticals': 'pharma',
        'drugs':           'pharma',
        'pharmacological': 'pharma',
        'pharmacotherapy': 'pharma',
        'treatment':       'pharma',
        'treatments':      'pharma',
        # subtypes plural/synonyms
        'subtype':         'subtypes',
        'endotype':        'subtypes',
        'endotypes':       'subtypes',
        'heterogeneity':   'subtypes',
        # mechanisms/pathways singular
        'mechanism':       'mechanisms',
        'pathway':         'pathways',
        # protective synonyms
        'lifestyle':       'protective',
        'prevention':      'protective',
        'protective_factors': 'protective',
        # comorbidities typos only
        'comorbidity':     'comorbidities',
        'co_morbidities':  'comorbidities',
        'comorbities':     'comorbidities',   # common Zotero typo
    }
    return _ALIASES.get(slug, slug)


# ── PubMed E-utilities DOI→PMID resolver ──────────────────────────────────────

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def _doi_to_pmid_bulk(dois: list[str],
                      ncbi_api_key: str | None = None) -> dict[str, str]:
    """
    Resolve a list of DOIs to PMIDs using PubMed E-utilities esearch.

    Makes one request per DOI (the batch-search approach doesn't return
    per-DOI mappings).  Rate: ~3 req/s without NCBI key, ~10/s with key.

    Parameters
    ----------
    dois          : list of DOI strings (may include 'doi:' prefix)
    ncbi_api_key  : optional NCBI API key for higher rate limits

    Returns
    -------
    dict  {doi_lowercase: pmid_string}  — only entries where a PMID was found
    """
    result: dict[str, str] = {}
    if not dois:
        return result

    print(f"  Resolving {len(dois)} DOIs → PMIDs via PubMed E-utilities ...")
    for doi in dois:
        doi_clean = re.sub(r'^doi:\s*', '', doi.strip(), flags=re.IGNORECASE)
        if not doi_clean:
            continue
        params: dict = {
            'db': 'pubmed',
            'term': f'{doi_clean}[doi]',
            'retmode': 'json',
            'retmax': 1,
        }
        if ncbi_api_key:
            params['api_key'] = ncbi_api_key

        data = _get(f"{EUTILS_BASE}/esearch.fcgi", params=params)
        time.sleep(0.35)   # stay within ~3 req/s without key

        if data:
            ids = data.get('esearchresult', {}).get('idlist', [])
            if ids:
                result[doi_clean.lower()] = ids[0]

    found = len(result)
    print(f"  DOI→PMID: resolved {found}/{len(dois)}")
    return result


# ── PubTator3 entity extraction ────────────────────────────────────────────────

def _fetch_pubtator3_annotations(pmids: list[str]) -> list[dict]:
    """
    Fetch entity annotations from PubTator3 biocjson endpoint in batches.
    Returns a flat list of {'text', 'type', 'pmid'} dicts.
    """
    entities: list[dict] = []

    for i in range(0, len(pmids), PUBTATOR3_BATCH_SIZE):
        batch = pmids[i:i + PUBTATOR3_BATCH_SIZE]
        pmid_str = ','.join(batch)
        data = _get(
            f"{PUBTATOR3_BASE}/publications/export/biocjson",
            params={'pmids': pmid_str, 'full': 'true'},
        )
        time.sleep(0.4)
        if not data:
            continue

        # BioC JSON: data['PubTator3'] is a list of publication dicts
        for pub in data.get('PubTator3', []):
            pmid = str(pub.get('pmid', ''))
            for passage in pub.get('passages', []):
                for ann in passage.get('annotations', []):
                    infons = ann.get('infons', {})
                    etype  = infons.get('type', '').strip()
                    etext  = ann.get('text', '').strip()
                    if etext and etype and etype not in EXCLUDED_ENTITY_TYPES:
                        entities.append({'text': etext, 'type': etype, 'pmid': pmid})

        print(f"    PubTator3: {len(batch)} papers in batch "
              f"{i // PUBTATOR3_BATCH_SIZE + 1} → "
              f"{sum(1 for e in entities if e['pmid'] in batch)} entities")

    return entities


def _aggregate_entities(entities: list[dict]) -> list[dict]:
    """
    Deduplicate and count entities by (text, type). Returns list sorted by freq.
    """
    counts:  Counter     = Counter()
    etypes:  dict        = {}
    for e in entities:
        key = e['text'].strip()
        if key:
            counts[key] += 1
            etypes.setdefault(key, e['type'])

    return [
        {'text': t, 'type': etypes[t], 'freq': c}
        for t, c in counts.most_common()
        if c >= MIN_ENTITY_FREQ
    ]


# ── Claude API classification ──────────────────────────────────────────────────

def _get_category_description(cat: str, phenotype: str) -> str:
    """
    Return a human-readable description for a category slug.
    Uses CATEGORY_DESCRIPTIONS when available; generates a generic fallback
    for dynamically discovered categories not in the default list.
    """
    if cat in CATEGORY_DESCRIPTIONS:
        return CATEGORY_DESCRIPTIONS[cat]
    readable = cat.replace('_', ' ')
    return f"terms related to {readable} aspects of {phenotype}"


def _build_classification_prompt(phenotype: str, n_papers: int,
                                  min_freq: int, entity_list: str,
                                  categories: list[str] | None = None,
                                  category_hint: str | None = None) -> str:
    """
    Build the Claude entity-classification prompt.

    Uses the provided `categories` list (dynamically derived from Zotero
    sub-folders) so the prompt always reflects the actual folder structure
    rather than a hardcoded schema.  Falls back to SEED_CATEGORIES when
    `categories` is None (flat collection / PMID file workflow).
    """
    cats = categories or SEED_CATEGORIES
    cat_lines = '\n'.join(
        f"- {cat:<18}: {_get_category_description(cat, phenotype)}"
        for cat in cats
    )
    example_json = json.dumps(
        {cat: ['example term'] for cat in cats},
        indent=2,
    )

    hint_note = (
        f"\n\nContext: these papers specifically discuss the '{category_hint}' "
        f"aspect of {phenotype} — weight entity assignments accordingly. "
        f"Terms that clearly belong to '{category_hint}' should be placed there."
        if category_hint else ""
    )

    return f"""\
You are a biomedical expert building a structured knowledge base about {phenotype}.

The entities below were extracted from {n_papers} research papers using \
PubTator3 entity recognition. Each entity appears in at least {min_freq} papers.
Classify each entity into the most appropriate category for the knowledge base.

Categories:
{cat_lines}
- exclude             : entity is not relevant to {phenotype}, too generic, a \
species/cell-line name, or a common English word

Extracted entities (entity text | entity type | frequency):
{entity_list}

Rules:
1. Return ONLY a JSON object — no preamble, no markdown fences, no comments.
2. Keys are exactly the {len(cats)} category names listed above (NOT 'exclude').
3. Values are arrays of standardised term strings.
4. Standardise terminology: use English, spell out Greek letters \
(e.g. "beta-cell dysfunction" not "β-cell dysfunction").
5. Omit terms that are irrelevant to {phenotype} or too generic.
6. Each term may appear in AT MOST ONE category.
7. Keep granular, specific terms (e.g. "SIDD", "HbA1c", "DRB5*0301") \
rather than collapsing to generic labels.{hint_note}

Return format (example):
{example_json}
"""


def classify_with_claude(
    entities: list[dict],
    phenotype: str,
    n_papers: int,
    api_key: str | None = None,
    category_hint: str | None = None,
    categories: list[str] | None = None,
) -> dict[str, list[str]]:
    """
    Use Claude claude-opus-4-6 with adaptive thinking to classify extracted entities
    into seed categories.

    Parameters
    ----------
    entities      : aggregated entity list from _aggregate_entities()
    phenotype     : primary phenotype string
    n_papers      : total number of source papers (for prompt context)
    api_key       : Anthropic API key (falls back to ANTHROPIC_API_KEY env var)
    category_hint : if not None, weight classification towards this category
    categories    : category list derived from Zotero sub-folders;
                    falls back to SEED_CATEGORIES when None

    Returns
    -------
    dict  {category: [term, ...]}  — may be empty if classification fails
    """
    if not _ANTHROPIC_AVAILABLE:
        print("  [WARN] 'anthropic' package not installed. "
              "Run: pip install anthropic\n"
              "  Returning raw entity list without classification.")
        return {}

    # Format entity list for prompt (top 250 max)
    top = entities[:250]
    entity_lines = '\n'.join(
        f"  {e['text']} | {e['type']} | freq={e['freq']}"
        for e in top
    )

    cats = categories or SEED_CATEGORIES

    prompt = _build_classification_prompt(
        phenotype=phenotype,
        n_papers=n_papers,
        min_freq=MIN_ENTITY_FREQ,
        entity_list=entity_lines,
        categories=cats,
        category_hint=category_hint,
    )

    client = _anthropic.Anthropic(api_key=api_key)

    print(f"  Sending {len(top)} entities to Claude for classification...")
    try:
        response = client.messages.create(
            model="claude-opus-4-6",
            max_tokens=16000,
            thinking={"type": "adaptive"},
            messages=[{"role": "user", "content": prompt}],
        )
    except Exception as exc:
        print(f"  [WARN] Claude API error: {exc}")
        return {}

    def _extract_json_text(content_blocks) -> str:
        """
        Pull the JSON string out of the response.

        Priority:
        1. Non-empty text blocks (normal case).
        2. Thinking block content — Claude sometimes places the final
           JSON answer inside the thinking block when the text block is
           empty (known edge-case with extended thinking).
        """
        # 1. Text blocks first
        for block in content_blocks:
            if block.type == "text" and block.text.strip():
                return block.text
        # 2. Thinking block fallback — find the last {...} in thinking content
        for block in content_blocks:
            if block.type == "thinking":
                raw = getattr(block, 'thinking', '') or ''
                # Grab the last JSON object in the thinking output
                matches = re.findall(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)?\}',
                                     raw, re.DOTALL)
                if matches:
                    return matches[-1]
        return ""

    text = _extract_json_text(response.content)

    # Strip any residual <thinking>…</thinking> wrappers and markdown fences
    text = re.sub(r'<thinking>.*?</thinking>\s*', '', text,
                  flags=re.DOTALL | re.IGNORECASE).strip()
    text = re.sub(r'^```[a-z]*\n?', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n?```$',        '', text, flags=re.MULTILINE)
    text = text.strip()

    try:
        result = json.loads(text)
        # Normalise: ensure all discovered categories are present
        return {cat: result.get(cat, []) for cat in cats}
    except json.JSONDecodeError as exc:
        print(f"  [WARN] Could not parse Claude response as JSON: {exc}")
        print(f"  Response (first 500 chars): {text[:500]}")
        return {}


# ── Seeds builder ──────────────────────────────────────────────────────────────

def _merge_seeds(existing: dict, new_terms: dict,
                 categories: list[str] | None = None) -> dict:
    """
    Merge new_terms into existing seeds dict.
    Existing curated terms are preserved; new terms are appended if not present.
    Uses the union of all categories seen across both dicts when categories
    is None, so merging never silently drops a category.
    """
    cats = categories or sorted(set(existing) | set(new_terms))
    merged = {cat: list(existing.get(cat, [])) for cat in cats}
    for cat in cats:
        existing_lower = {t.lower() for t in merged[cat]}
        for term in new_terms.get(cat, []):
            if term.lower() not in existing_lower:
                merged[cat].append(term)
                existing_lower.add(term.lower())
    return merged


def build_seeds_from_papers(
    papers_by_category: dict[str, list[dict]],
    phenotype: str,
    output_path: str,
    api_key: str | None = None,
    ncbi_api_key: str | None = None,
    merge: bool = False,
    categories: list[str] | None = None,
) -> dict:
    """
    Core routine: extract entities, classify, write seeds JSON.

    Parameters
    ----------
    papers_by_category : dict mapping category name → list of paper dicts.
                         Use 'all' as the key if papers are not pre-categorised.
    phenotype          : primary phenotype string
    output_path        : output JSON file path
    api_key            : Anthropic API key
    ncbi_api_key       : optional NCBI API key for higher E-utilities rate limits
    merge              : if True, merge with existing seeds file
    categories         : ordered category list derived from Zotero sub-folders;
                         falls back to SEED_CATEGORIES for flat/PMID workflows

    Returns
    -------
    dict  final seeds dict
    """
    all_papers: list[dict] = [p for papers in papers_by_category.values()
                               for p in papers]
    n_total = len(all_papers)
    print(f"\n  Total papers : {n_total}")

    # ── DOI→PMID fallback for papers missing a PMID ───────────────────────
    missing_pmid = [p for p in all_papers if not p['pmid'] and p.get('doi')]
    if missing_pmid:
        print(f"  [INFO] {len(missing_pmid)} papers missing PMID but have DOI — "
              f"attempting DOI→PMID lookup ...")
        doi_map = _doi_to_pmid_bulk(
            [p['doi'] for p in missing_pmid],
            ncbi_api_key=ncbi_api_key,
        )
        # Patch the PMID in every matching paper (across all categories)
        for cat_papers in papers_by_category.values():
            for p in cat_papers:
                if not p['pmid'] and p.get('doi'):
                    resolved = doi_map.get(p['doi'].lower().lstrip('doi:').strip())
                    if resolved:
                        p['pmid'] = resolved

    # ── PubTator3 entity extraction ────────────────────────────────────────
    all_pmids = [p['pmid'] for p in all_papers if p['pmid']]
    n_no_pmid = n_total - len(all_pmids)
    print(f"  With PMID    : {len(all_pmids)}  (without: {n_no_pmid})")

    if n_no_pmid:
        print(f"  [NOTE] {n_no_pmid} papers still have no PMID — "
              f"their abstracts won't be entity-annotated.")

    # Use caller-supplied category list (from Zotero folders) or the default
    cats = categories or SEED_CATEGORIES
    combined_seeds: dict[str, list[str]] = {cat: [] for cat in cats}

    if 'all' in papers_by_category or len(papers_by_category) == 1:
        # Flat collection: extract all entities then classify
        key   = list(papers_by_category.keys())[0]
        plist = papers_by_category[key]
        pmids = [p['pmid'] for p in plist if p['pmid']]
        print(f"\n[entity extraction]  ({len(pmids)} papers with PMID)")

        raw_entities = _fetch_pubtator3_annotations(pmids)
        aggregated   = _aggregate_entities(raw_entities)
        print(f"  Unique entities (freq≥{MIN_ENTITY_FREQ}): {len(aggregated)}")

        if aggregated:
            classified = classify_with_claude(
                aggregated, phenotype, len(plist), api_key,
                category_hint=key if key != 'all' else None,
                categories=cats,
            )
            combined_seeds = classified if classified else combined_seeds
        else:
            print("  [WARN] No entities extracted — output seeds will be empty.")

    else:
        # Sub-collection → category mapping: classify per category
        for category, plist in papers_by_category.items():
            pmids = [p['pmid'] for p in plist if p['pmid']]
            if not pmids:
                print(f"\n[{category}]  no papers with PMID — skipping")
                continue

            print(f"\n[{category}]  {len(pmids)} papers")
            raw_entities = _fetch_pubtator3_annotations(pmids)
            aggregated   = _aggregate_entities(raw_entities)
            print(f"  Unique entities (freq≥{MIN_ENTITY_FREQ}): {len(aggregated)}")

            if aggregated:
                classified = classify_with_claude(
                    aggregated, phenotype, len(plist), api_key,
                    category_hint=category,
                    categories=cats,
                )
                # Prefer the matched category; fall back to any non-empty result
                if classified.get(category):
                    combined_seeds[category].extend(classified[category])
                else:
                    for cat, terms in classified.items():
                        if cat in combined_seeds:
                            combined_seeds[cat].extend(terms)

    # Deduplicate within each category
    for cat in cats:
        seen:  set  = set()
        dedup: list = []
        for t in combined_seeds[cat]:
            tl = t.lower()
            if tl not in seen:
                dedup.append(t)
                seen.add(tl)
        combined_seeds[cat] = dedup

    # ── Merge with existing file if requested ──────────────────────────────
    if merge and os.path.exists(output_path):
        try:
            with open(output_path) as f:
                existing = json.load(f)
            existing_seeds = {
                k: v for k, v in existing.items()
                if isinstance(v, list)
            }
            combined_seeds = _merge_seeds(existing_seeds, combined_seeds,
                                          categories=cats)
            print(f"\n  Merged with existing seeds at {output_path}")
        except Exception as exc:
            print(f"  [WARN] Could not read existing seeds for merge: {exc}")

    # ── Write output ───────────────────────────────────────────────────────
    output = {
        'phenotype':  phenotype,
        '_generated': datetime.now().isoformat(),
        '_note':      (
            "Auto-generated from literature. Review and curate before use. "
            "Add any missing terms from disease-specific knowledge."
        ),
    }
    output.update(combined_seeds)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    total_terms = sum(len(v) for v in combined_seeds.values())
    print(f"\n  Seeds written: {output_path}")
    print(f"  Total terms  : {total_terms}")
    for cat in cats:
        n = len(combined_seeds.get(cat, []))
        sample = combined_seeds.get(cat, [])[:4]
        print(f"    {cat:<18} {n:>3} terms   {sample}")

    return combined_seeds


# ── Zotero workflow ────────────────────────────────────────────────────────────

def run_from_zotero(
    phenotype: str,
    collection_path: str,
    output_path: str,
    db_path: str = DEFAULT_ZOTERO_DB,
    category_map: dict[str, str] | None = None,
    api_key: str | None = None,
    ncbi_api_key: str | None = None,
    merge: bool = False,
) -> dict:
    """
    Build seeds from a Zotero collection.

    Parameters
    ----------
    phenotype       : primary phenotype string
    collection_path : Zotero collection path.
                      Use slash notation to disambiguate collections with
                      the same name at different levels:
                        'seedCollection/T2D'  — unambiguous (recommended)
                        'T2D'                 — matches any collection named T2D
    output_path     : output JSON file path
    db_path         : path to zotero.sqlite
    category_map    : explicit {subcollection_name: category} mapping;
                      if None, auto-detect from subcollection names via
                      _folder_name_to_category_slug()
    api_key         : Anthropic API key
    ncbi_api_key    : optional NCBI API key for faster DOI→PMID resolution
    merge           : merge with existing seeds file
    """
    print(f"\nZotero database  : {db_path}")
    print(f"Seed collection  : {collection_path}")
    if not os.path.exists(db_path):
        raise FileNotFoundError(
            f"Zotero database not found at '{db_path}'. "
            f"Pass --zotero_db with the correct path."
        )

    conn = _open_zotero_db(db_path)
    try:
        # ── Resolve collection by path ─────────────────────────────────────
        col_ids = _resolve_collection_path(conn, collection_path)
        if not col_ids:
            raise ValueError(
                f"Collection '{collection_path}' not found in Zotero.\n"
                f"Run --list_collections to see available collections and "
                f"their full paths."
            )
        if len(col_ids) > 1:
            print(f"  [NOTE] '{collection_path}' matched {len(col_ids)} collections "
                  f"— papers from all will be included.  Use a slash path "
                  f"(e.g. 'seedCollection/T2D') to match exactly one.")

        # ── Check for sub-collections ─────────────────────────────────────
        subcols = _get_subcollections_by_ids(conn, col_ids)
        papers_by_category: dict[str, list[dict]] = {}

        if subcols:
            print(f"\nSub-collections found under '{collection_path}': {len(subcols)}")
            for sc in subcols:
                sc_name = sc['collectionName']
                sc_id   = sc['collectionID']
                # Include the sub-collection and any of its own children
                sc_all_ids = _expand_to_descendants(conn, [sc_id])
                papers     = _query_papers_in_collections(conn, sc_all_ids)

                # Category: explicit override → slug derived from folder name
                # Any slug is valid — categories are defined by the folders,
                # not by a hardcoded list.
                if category_map and sc_name in category_map:
                    category = category_map[sc_name]
                else:
                    category = _folder_name_to_category_slug(sc_name)

                papers_by_category.setdefault(category, []).extend(papers)
                print(f"  '{sc_name}' → '{category}'  ({len(papers)} papers)")

        # Also include papers directly in the parent collection(s)
        parent_papers = _query_papers_in_collections(conn, col_ids)
        if parent_papers:
            papers_by_category.setdefault('all', []).extend(parent_papers)
            print(f"\nDirect papers in '{collection_path}': {len(parent_papers)}")

        if not papers_by_category:
            raise ValueError(
                f"No papers found in '{collection_path}'.  "
                f"Check the collection path with --list_collections."
            )

    finally:
        conn.close()

    # Preserve the order in which sub-folders appeared; exclude the 'all' key
    # which is used for papers sitting directly in the parent collection.
    discovered_categories = [c for c in papers_by_category if c != 'all'] or None

    return build_seeds_from_papers(
        papers_by_category, phenotype, output_path,
        api_key=api_key, ncbi_api_key=ncbi_api_key, merge=merge,
        categories=discovered_categories,
    )


# ── PMID file workflow ─────────────────────────────────────────────────────────

def run_from_pmids(
    phenotype: str,
    pmids_file: str,
    output_path: str,
    api_key: str | None = None,
    ncbi_api_key: str | None = None,
    merge: bool = False,
) -> dict:
    """
    Build seeds from a plain-text file of PMIDs (one per line).
    """
    with open(pmids_file) as f:
        pmids = [line.strip() for line in f if line.strip().isdigit()]

    if not pmids:
        raise ValueError(f"No valid PMIDs found in '{pmids_file}'.")

    print(f"  PMIDs loaded: {len(pmids)}")
    papers = [{'title': '', 'abstract': '', 'doi': '', 'pmid': p} for p in pmids]
    return build_seeds_from_papers(
        {'all': papers}, phenotype, output_path,
        api_key=api_key, ncbi_api_key=ncbi_api_key, merge=merge,
    )


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            "Build a phenotype seeds JSON from Zotero collections or a PMID list.\n\n"
            "Extracts entity annotations via PubTator3 and classifies them into\n"
            "seed categories using Claude claude-opus-4-6 (adaptive thinking).\n\n"
            "Output is written to a seeds JSON file usable by\n"
            "phenotype_background_analysis.py."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        '--phenotype', default=None,
        help="Primary phenotype (e.g. 'type 2 diabetes')",
    )
    parser.add_argument(
        '--output', default=None,
        help=(
            "Output JSON file path.\n"
            "If omitted, defaults to $PHENO_PATH/phenotype_seeds.json\n"
            "(following the project PHENO_PATH convention).\n"
            "Example: scripts/phenotype_seeds/type_2_diabetes.json"
        ),
    )
    parser.add_argument(
        '--zotero_collection', default=None,
        help=(
            "Zotero collection path for the SEED articles.\n"
            "Use slash notation to pinpoint the right collection when\n"
            "multiple collections share the same name:\n"
            "  'seedCollection/T2D'       — unambiguous (recommended)\n"
            "  'T2D'                      — matches any collection named T2D\n"
            "The named collection's sub-folders define seed categories;\n"
            "sub-folder names are mapped to category slugs automatically."
        ),
    )
    parser.add_argument(
        '--zotero_db', default=DEFAULT_ZOTERO_DB,
        help=f"Path to Zotero SQLite database (default: {DEFAULT_ZOTERO_DB})",
    )
    parser.add_argument(
        '--category_map', default=None,
        help=(
            "Explicit sub-collection → category mapping, comma-separated.\n"
            "Example: \"T2D Subtypes:subtypes,T2D Pathways:pathways\""
        ),
    )
    parser.add_argument(
        '--pmids_file', default=None,
        help="Path to a text file with one PMID per line (alternative to Zotero)",
    )
    parser.add_argument(
        '--api_key', default=None,
        help=(
            "Anthropic API key for Claude classification.\n"
            "Falls back to ANTHROPIC_API_KEY environment variable.\n"
            "If neither is set, raw entity lists are written without classification."
        ),
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help=(
            "NCBI API key for faster DOI→PMID resolution via E-utilities\n"
            "(raises rate limit from ~3 to ~10 req/s).\n"
            "Get a free key at: https://www.ncbi.nlm.nih.gov/account/"
        ),
    )
    parser.add_argument(
        '--merge', action='store_true',
        help="Merge new terms into existing seeds file rather than overwrite",
    )
    parser.add_argument(
        '--list_collections', action='store_true',
        help="List all Zotero collections and exit (no other args needed)",
    )

    args = parser.parse_args()

    # ── List mode ─────────────────────────────────────────────────────────
    if args.list_collections:
        db = args.zotero_db
        if not os.path.exists(db):
            print(f"Zotero database not found at '{db}'.")
            print(f"Pass --zotero_db with the correct path.")
            raise SystemExit(1)

        cols = list_collections(db)
        print(f"\nZotero collections in {db}:\n")
        print(f"  {'FULL PATH':<55}  ITEMS")
        print(f"  {'-'*55}  -----")
        for c in cols:
            print(f"  {c['full_path']:<55}  {c['item_count']}")
        print(f"\n  Use the FULL PATH with --zotero_collection, e.g.:")
        print(f"    --zotero_collection 'seedCollection/T2D'")
        print(f"  Use the root-level name with --zotero_corpus, e.g.:")
        print(f"    --zotero_corpus T2D")
        raise SystemExit(0)

    # ── Validate args ─────────────────────────────────────────────────────
    if not args.phenotype:
        parser.error("--phenotype is required")
    if not args.zotero_collection and not args.pmids_file:
        parser.error("Either --zotero_collection or --pmids_file is required")

    # Resolve output path: CLI → PHENO_PATH env var → error
    output_path = args.output
    if not output_path:
        pheno_path_env = os.environ.get("PHENO_PATH")
        if pheno_path_env:
            output_path = os.path.join(pheno_path_env, "phenotype_seeds.json")
        else:
            parser.error(
                "--output is required (or set the PHENO_PATH environment variable, "
                "e.g. export PHENO_PATH=$PROJECT_ROOT/results/type_2_diabetes)"
            )

    print("=" * 60)
    print("BUILD SEEDS FROM ARTICLES")
    print("=" * 60)
    print(f"  Phenotype : {args.phenotype}")
    print(f"  Output    : {output_path}")

    # Parse category map
    cat_map: dict[str, str] | None = None
    if args.category_map:
        cat_map = {}
        for pair in args.category_map.split(','):
            pair = pair.strip()
            if ':' in pair:
                sc_name, category = pair.split(':', 1)
                cat_map[sc_name.strip()] = category.strip()

    if args.zotero_collection:
        run_from_zotero(
            phenotype=args.phenotype,
            collection_path=args.zotero_collection,
            output_path=output_path,
            db_path=args.zotero_db,
            category_map=cat_map,
            api_key=args.api_key,
            ncbi_api_key=args.ncbi_api_key,
            merge=args.merge,
        )
    else:
        run_from_pmids(
            phenotype=args.phenotype,
            pmids_file=args.pmids_file,
            output_path=output_path,
            api_key=args.api_key,
            ncbi_api_key=args.ncbi_api_key,
            merge=args.merge,
        )
