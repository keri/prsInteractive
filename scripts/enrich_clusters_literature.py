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
  - Zotero corpus       — optional; full-text PDFs + abstracts from a curated
                          local Zotero collection; Claude claude-opus-4-6 sentence-level
                          sentiment analysis (RISK / PROTECTIVE / NEUTRAL)

Query strategy per gene/SNP × phenotype
----------------------------------------
1. DisGeNET GDA score (gene-disease association)
2. Open Targets association score (integrated evidence)
3. PubTator3 entity-annotated co-occurrence count
4. PubMed keyword co-occurrence count (fallback / complement)
5. ClinVar pathogenicity (SNPs and short gene symbols only)
6. Zotero corpus sentiment score (optional; 20% weight when active)
7. Combined score = weighted sum of the above

Usage
-----
  python enrich_clusters_literature.py \\
    --cluster_dir /path/to/cohortBnmf/cardio \\
    --phenotypes "type 2 diabetes,epilepsy,cardio" \\
    --ncbi_api_key YOUR_NCBI_KEY \\
    [--disgenet_api_key YOUR_DISGENET_KEY] \\
    [--zotero_corpus T2D --anthropic_api_key YOUR_ANTHROPIC_KEY]
"""

import argparse
import json
import math
import os
import re
import sqlite3
import time
from collections import defaultdict

import pandas as pd
import requests

# ── Optional: Anthropic SDK for Zotero corpus sentiment analysis ───────────────
try:
    import anthropic as _anthropic
    _ANTHROPIC_AVAILABLE = True
except ImportError:
    _ANTHROPIC_AVAILABLE = False

# ── Optional: PDF text extraction ─────────────────────────────────────────────
try:
    import pdfplumber as _pdfplumber
    _PDFPLUMBER_AVAILABLE = True
except ImportError:
    _PDFPLUMBER_AVAILABLE = False

try:
    import pypdf as _pypdf
    _PYPDF_AVAILABLE = True
except ImportError:
    _PYPDF_AVAILABLE = False


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

# Default Zotero storage directory
DEFAULT_ZOTERO_DB      = os.path.expanduser('~/Zotero/zotero.sqlite')
DEFAULT_ZOTERO_STORAGE = os.path.expanduser('~/Zotero/storage')

# Max sentences sent to Claude per (feature, phenotype) pair — controls cost
CORPUS_MAX_SENTENCES = 120

# Weight of Zotero corpus evidence in the combined score (when corpus is active)
WEIGHT_CORPUS = 0.20

# Rescaled weights when corpus is included (original total = 1.0)
# Corpus takes 0.20; remaining 0.80 distributed proportionally
_W_BASE_TOTAL = (
    WEIGHT_DISGENET + WEIGHT_OPEN_TARGETS + WEIGHT_PUBTATOR3 + WEIGHT_CLINVAR
)
WEIGHT_DISGENET_CORP     = WEIGHT_DISGENET     / _W_BASE_TOTAL * (1.0 - WEIGHT_CORPUS)
WEIGHT_OPEN_TARGETS_CORP = WEIGHT_OPEN_TARGETS / _W_BASE_TOTAL * (1.0 - WEIGHT_CORPUS)
WEIGHT_PUBTATOR3_CORP    = WEIGHT_PUBTATOR3    / _W_BASE_TOTAL * (1.0 - WEIGHT_CORPUS)
WEIGHT_CLINVAR_CORP      = WEIGHT_CLINVAR      / _W_BASE_TOTAL * (1.0 - WEIGHT_CORPUS)


# ── Zotero corpus loader ───────────────────────────────────────────────────────

class ZoteroCorpus:
    """
    Local Zotero literature corpus for feature-level enrichment.

    Loads all papers from a named Zotero collection, extracts full text from
    PDF attachments where available (via pdfplumber → pypdf fallback), and
    falls back to PubTator3-fetched abstracts for papers without local PDFs.

    Usage
    -----
    corpus = ZoteroCorpus('T2D', db_path='~/Zotero/zotero.sqlite')
    corpus.load()                                     # one-time load
    matches = corpus.search(['KAZN', 'kazrin'])       # list of matching papers
    """

    def __init__(self,
                 collection_name: str,
                 db_path: str = DEFAULT_ZOTERO_DB,
                 storage_dir: str = DEFAULT_ZOTERO_STORAGE):
        self.collection_name = collection_name
        self.db_path         = os.path.expanduser(db_path)
        self.storage_dir     = os.path.expanduser(storage_dir)
        self.papers: dict[str, dict] = {}   # corpus_key → paper dict
        self._loaded = False

    # ── SQLite helpers ─────────────────────────────────────────────────────

    def _open_db(self) -> sqlite3.Connection:
        uri = f'file:{self.db_path}?mode=ro&immutable=1'
        conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        return conn

    def _get_collection_ids(self, conn: sqlite3.Connection) -> list[int]:
        """Recursive CTE — returns the named collection + all descendants."""
        root_rows = conn.execute(
            "SELECT collectionID FROM collections "
            "WHERE LOWER(collectionName) = LOWER(?)",
            (self.collection_name,),
        ).fetchall()
        if not root_rows:
            root_rows = conn.execute(
                "SELECT collectionID FROM collections "
                "WHERE LOWER(collectionName) LIKE LOWER(?)",
                (f'%{self.collection_name}%',),
            ).fetchall()
        if not root_rows:
            return []
        root_ids = [r[0] for r in root_rows]
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

    def _query_papers(self, conn: sqlite3.Connection,
                      collection_ids: list[int]) -> list[dict]:
        if not collection_ids:
            return []
        placeholders = ','.join('?' * len(collection_ids))
        rows = conn.execute(f"""
            SELECT DISTINCT
                i.itemID,
                i.key                                                          AS item_key,
                it.typeName                                                    AS item_type,
                MAX(CASE WHEN f.fieldName = 'title'        THEN idv.value END) AS title,
                MAX(CASE WHEN f.fieldName = 'abstractNote' THEN idv.value END) AS abstract,
                MAX(CASE WHEN f.fieldName = 'DOI'          THEN idv.value END) AS doi,
                MAX(CASE WHEN f.fieldName = 'PMID'         THEN idv.value END) AS pmid_field,
                MAX(CASE WHEN f.fieldName = 'PMCID'        THEN idv.value END) AS pmcid_field,
                MAX(CASE WHEN f.fieldName = 'extra'        THEN idv.value END) AS extra,
                MAX(CASE WHEN f.fieldName = 'url'          THEN idv.value END) AS url
            FROM items i
            JOIN itemTypes it   ON i.itemTypeID  = it.itemTypeID
            JOIN collectionItems ci ON i.itemID  = ci.itemID
            LEFT JOIN itemData id    ON i.itemID  = id.itemID
            LEFT JOIN fields f       ON id.fieldID = f.fieldID
            LEFT JOIN itemDataValues idv ON id.valueID = idv.valueID
            WHERE ci.collectionID IN ({placeholders})
              AND it.typeName NOT IN ('note', 'attachment')
            GROUP BY i.itemID
        """, collection_ids).fetchall()
        return [dict(r) for r in rows]

    def _get_pdf_paths(self, conn: sqlite3.Connection,
                       item_ids: list[int]) -> dict[int, str]:
        """Return {itemID: local_pdf_path} for items that have a stored PDF."""
        if not item_ids:
            return {}
        placeholders = ','.join('?' * len(item_ids))
        rows = conn.execute(f"""
            SELECT
                ia.parentItemID,
                att.key       AS att_key,
                ia.path       AS attachment_path,
                ia.contentType
            FROM itemAttachments ia
            JOIN items att ON ia.itemID = att.itemID
            WHERE ia.parentItemID IN ({placeholders})
              AND ia.contentType = 'application/pdf'
            ORDER BY ia.itemID
        """, item_ids).fetchall()

        pdf_map: dict[int, str] = {}
        for row in rows:
            parent_id = row['parentItemID']
            if parent_id in pdf_map:
                continue   # already have a PDF for this item
            att_key = row['att_key']
            raw_path = row['attachment_path'] or ''
            if raw_path.startswith('storage:'):
                filename = raw_path[8:]     # strip 'storage:' prefix
                pdf_path = os.path.join(self.storage_dir, att_key, filename)
            elif os.path.isabs(raw_path):
                pdf_path = raw_path         # linked file: absolute path
            else:
                # Fallback: att_key directory with whatever filename we have
                pdf_path = os.path.join(self.storage_dir, att_key, raw_path)
            pdf_map[parent_id] = pdf_path
        return pdf_map

    # ── PDF text extraction ────────────────────────────────────────────────

    @staticmethod
    def _extract_pdf_text(pdf_path: str, max_pages: int = 60) -> str:
        """Extract plain text from a PDF. Returns '' if extraction fails."""
        if not pdf_path or not os.path.exists(pdf_path):
            return ''

        # Preferred: pdfplumber (handles complex layouts well)
        if _PDFPLUMBER_AVAILABLE:
            try:
                with _pdfplumber.open(pdf_path) as pdf:
                    pages = pdf.pages[:max_pages]
                    return ' '.join(
                        t for page in pages
                        for t in [page.extract_text() or '']
                        if t
                    )
            except Exception:
                pass

        # Fallback: pypdf
        if _PYPDF_AVAILABLE:
            try:
                with open(pdf_path, 'rb') as fh:
                    reader = _pypdf.PdfReader(fh)
                    return ' '.join(
                        reader.pages[i].extract_text() or ''
                        for i in range(min(max_pages, len(reader.pages)))
                    )
            except Exception:
                pass

        return ''

    # ── Corpus loading ─────────────────────────────────────────────────────

    @staticmethod
    def _extract_pmid(extra: str, url: str) -> str:
        if extra:
            m = re.search(r'PMID:\s*(\d{6,})', extra, re.IGNORECASE)
            if m:
                return m.group(1)
        if url:
            m = re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d{6,})', url)
            if m:
                return m.group(1)
        return ''

    @staticmethod
    def _extract_doi(doi_field: str, extra: str, url: str) -> str:
        """Extract DOI from Zotero fields (DOI field → extra → URL)."""
        if doi_field:
            return doi_field.strip()
        if extra:
            m = re.search(r'DOI:\s*(10\.\S+)', extra, re.IGNORECASE)
            if m:
                return m.group(1).rstrip('.,;)')
        if url and '10.' in url:
            m = re.search(r'(10\.\d{4,}/\S+)', url)
            if m:
                return m.group(1).rstrip('.,;)')
        return ''

    # ── PMC full-text fetch ────────────────────────────────────────────────

    @staticmethod
    def _fetch_pmc_fulltext(pmid: str, api_key: str | None = None,
                             cache: dict | None = None,
                             doi: str = '',
                             pmcid: str = '') -> str:
        """
        Fetch full text from PubMed Central given any combination of identifiers.

        Resolution order (first that yields a PMCID wins):
        A. pmcid supplied directly — skip lookup entirely (fastest)
        B. pmid supplied  → elink (PMID→PMCID) via E-utilities
        C. doi supplied   → idconv API (DOI→PMCID) — most reliable for papers
           added to Zotero via DOI rather than PubMed import
        Then:
        D. efetch (PMCID → full-text XML)
        E. Parse body text, strip XML tags

        Returns plain text or '' if the paper is not in PMC / not open access.
        Respects NCBI rate limits (10 req/s with key, 3 req/s without).
        """
        if not pmid and not doi and not pmcid:
            return ''

        cache_key = f'pmc_{pmcid or pmid or doi}'
        if cache is not None and cache_key in cache:
            return cache[cache_key]

        _headers = {'User-Agent': 'ZoteroCorpus/1.0 (research; contact via GitHub)'}
        _base    = {'api_key': api_key} if api_key else {}
        _sleep   = 0.12 if api_key else 0.4
        resolved_pmcid = pmcid   # may already be set

        # ── Path A: direct PMCID — no lookup needed ────────────────────────
        # (resolved_pmcid already set above)

        # ── Path B: PMID → PMCID via elink ────────────────────────────────
        if pmid and not resolved_pmcid:
            try:
                r = requests.get(
                    f'{NCBI_BASE}/elink.fcgi',
                    params={**_base, 'dbfrom': 'pubmed', 'db': 'pmc',
                            'id': pmid, 'retmode': 'json'},
                    headers=_headers, timeout=15,
                )
                time.sleep(_sleep)
                if r.status_code == 200:
                    data = r.json()
                    for ls in data.get('linksets', []):
                        for ld in ls.get('linksetdbs', []):
                            if ld.get('dbto') == 'pmc':
                                links = ld.get('links', [])
                                if links:
                                    resolved_pmcid = str(links[0])
                                    break
            except Exception:
                pass

        # ── Path C: DOI → PMCID via idconv ────────────────────────────────
        if doi and not resolved_pmcid:
            try:
                r = requests.get(
                    'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/',
                    params={'ids': doi, 'format': 'json',
                            'tool': 'ZoteroCorpus', 'email': 'research@github'},
                    headers=_headers, timeout=15,
                )
                time.sleep(_sleep)
                if r.status_code == 200:
                    data = r.json()
                    for rec in data.get('records', []):
                        if rec.get('pmcid'):
                            resolved_pmcid = rec['pmcid']  # e.g. 'PMC1234567'
                            break
            except Exception:
                pass

        if not resolved_pmcid:
            if cache is not None:
                cache[cache_key] = ''
            return ''

        # ── Fetch full-text XML ────────────────────────────────────────────
        try:
            r = requests.get(
                f'{NCBI_BASE}/efetch.fcgi',
                params={**_base, 'db': 'pmc', 'id': resolved_pmcid,
                        'rettype': 'full', 'retmode': 'xml'},
                headers=_headers, timeout=45,
            )
            time.sleep(_sleep)
            if r.status_code != 200:
                if cache is not None:
                    cache[cache_key] = ''
                return ''
            xml = r.text
        except Exception:
            return ''

        # ── Extract body text ──────────────────────────────────────────────
        body = re.search(r'<body[^>]*>(.*?)</body>', xml, re.DOTALL)
        raw  = body.group(1) if body else xml
        # Remove front/back matter if body tag absent
        raw  = re.sub(r'<(?:front|back)>.*?</(?:front|back)>', '', raw,
                      flags=re.DOTALL)
        text = re.sub(r'<[^>]+>', ' ', raw)
        text = re.sub(r'\s+', ' ', text).strip()

        if cache is not None:
            cache[cache_key] = text
        return text

    # ── Unpaywall full-text fetch ──────────────────────────────────────────

    @staticmethod
    def _fetch_unpaywall_pdf(doi: str, email: str,
                              cache: dict | None = None) -> str:
        """
        Query Unpaywall for a legal open-access PDF URL, download and extract.

        Uses api.unpaywall.org/v2/{doi}?email={email}.
        Prefers url_for_pdf from best_oa_location; falls back to other OA
        locations.  Downloads PDF to a temp file, extracts text, then deletes.

        Returns plain text or '' if no open-access version is available.
        """
        if not doi or not email:
            return ''

        cache_key = f'unpaywall_{doi}'
        if cache is not None and cache_key in cache:
            return cache[cache_key]

        # ── Query Unpaywall ────────────────────────────────────────────────
        try:
            r = requests.get(
                f'https://api.unpaywall.org/v2/{doi}',
                params={'email': email},
                headers={'User-Agent': 'ZoteroCorpus/1.0 (research)'},
                timeout=15,
            )
            time.sleep(0.5)
            if r.status_code != 200:
                return ''
            data = r.json()
        except Exception:
            return ''

        # Find best open-access PDF URL
        pdf_url: str = ''
        best = data.get('best_oa_location') or {}
        pdf_url = best.get('url_for_pdf') or ''
        if not pdf_url:
            for loc in data.get('oa_locations', []):
                if loc.get('url_for_pdf'):
                    pdf_url = loc['url_for_pdf']
                    break
        if not pdf_url:
            if cache is not None:
                cache[cache_key] = ''
            return ''

        # ── Download PDF and extract text ──────────────────────────────────
        import tempfile
        try:
            r = requests.get(
                pdf_url, timeout=45,
                headers={'User-Agent': 'ZoteroCorpus/1.0 (research)'},
                allow_redirects=True,
            )
            ctype = r.headers.get('content-type', '').lower()
            if r.status_code != 200 or 'pdf' not in ctype:
                if cache is not None:
                    cache[cache_key] = ''
                return ''
            with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp:
                tmp.write(r.content)
                tmp_path = tmp.name
            text = ZoteroCorpus._extract_pdf_text(tmp_path)
            os.unlink(tmp_path)
        except Exception:
            text = ''

        if cache is not None:
            cache[cache_key] = text
        return text

    def load(self, ncbi_api_key: str | None = None,
             unpaywall_email: str | None = None) -> 'ZoteroCorpus':
        """
        Load all papers + PDF text from the Zotero collection.

        Call once before calling search(). Loading 200-500 papers with PDFs
        typically takes 30–120 seconds depending on PDF size and disk speed.

        For papers that only have abstracts (no local PDF), the method will
        attempt to fetch full text from:
          1. PubMed Central (requires ncbi_api_key for 10 req/s; 3/s without)
          2. Unpaywall open-access PDF (requires unpaywall_email)

        Results are cached in memory to avoid duplicate fetches within a run.
        """
        if self._loaded:
            return self

        if not os.path.exists(self.db_path):
            print(f"  [WARN] Zotero DB not found at '{self.db_path}' — "
                  f"corpus enrichment disabled.")
            self._loaded = True
            return self

        print(f"\n  Loading Zotero corpus: '{self.collection_name}' ...")
        conn = self._open_db()
        try:
            coll_ids = self._get_collection_ids(conn)
            if not coll_ids:
                print(f"  [WARN] Collection '{self.collection_name}' not found "
                      f"in Zotero — corpus enrichment disabled.")
                self._loaded = True
                return self

            papers = self._query_papers(conn, coll_ids)
            if not papers:
                print(f"  [WARN] No papers found in '{self.collection_name}'.")
                self._loaded = True
                return self

            item_ids  = [p['itemID'] for p in papers]
            pdf_map   = self._get_pdf_paths(conn, item_ids)
        finally:
            conn.close()

        # Build corpus entries, extract PDF text where available
        n_pdfs      = 0
        abs_only: list[str] = []   # keys of papers that ended up abstract-only

        for p in papers:
            # Prefer dedicated Zotero fields; fall back to parsing extra/url
            pmid     = (p.get('pmid_field') or '').strip() or \
                       self._extract_pmid(p.get('extra') or '', p.get('url') or '')
            pmcid    = (p.get('pmcid_field') or '').strip()
            # Normalise PMCID: some Zotero entries store 'PMC1234567', efetch
            # wants just the numeric part (or the full PMC-prefixed form is fine
            # for idconv; efetch prefers numeric, but PMC-prefixed also works).
            # Keep as-is — efetch accepts both 'PMC1234567' and '1234567'.
            doi      = self._extract_doi(p.get('doi') or '',
                                         p.get('extra') or '',
                                         p.get('url') or '')
            key      = pmid or f"item_{p['itemID']}"
            pdf_path = pdf_map.get(p['itemID'], '')
            full_text = ''
            is_full   = False
            if pdf_path:
                full_text = self._extract_pdf_text(pdf_path)
                if full_text:
                    is_full = True
                    n_pdfs += 1

            self.papers[key] = {
                'itemID':       p['itemID'],
                'pmid':         pmid,
                'pmcid':        pmcid,
                'doi':          doi,
                'title':        p.get('title') or '',
                'abstract':     p.get('abstract') or '',
                'full_text':    full_text,
                'is_full_text': is_full,
            }
            if not is_full:
                abs_only.append(key)

        # ── Fetch full text for abstract-only papers ───────────────────────
        # Shared in-memory cache so duplicate PMIDs/DOIs are not re-fetched.
        fetch_cache: dict = {}
        n_pmc = 0
        n_uw  = 0

        # Diagnostic: show identifier coverage before attempting fetches
        n_with_pmcid = sum(1 for k in abs_only if self.papers[k]['pmcid'])
        n_with_pmid  = sum(1 for k in abs_only if self.papers[k]['pmid'])
        n_with_doi   = sum(1 for k in abs_only if self.papers[k]['doi'])
        n_no_id      = sum(1 for k in abs_only
                           if not self.papers[k]['pmcid']
                           and not self.papers[k]['pmid']
                           and not self.papers[k]['doi'])
        print(f"    Abstract-only identifier coverage: "
              f"{n_with_pmcid} PMCID (direct fetch), "
              f"{n_with_pmid} PMID, {n_with_doi} DOI, "
              f"{n_no_id} none (title-only, skipped)")

        if (ncbi_api_key or unpaywall_email) and abs_only:
            try_pmc = bool(ncbi_api_key)
            try_uw  = bool(unpaywall_email)

            n_abs = len(abs_only)
            print(f"    Fetching full text for {n_abs} abstract-only papers "
                  f"(PMC={'yes (PMID+DOI→PMCID)' if try_pmc else 'no'}, "
                  f"Unpaywall={'yes' if try_uw else 'no'}) ...")

            for i, key in enumerate(abs_only, 1):
                paper = self.papers[key]
                text  = ''

                # ── PMC first: PMCID direct → PMID elink → DOI idconv ────
                if try_pmc and (paper['pmcid'] or paper['pmid'] or paper['doi']) \
                        and not text:
                    text = self._fetch_pmc_fulltext(
                        pmid=paper['pmid'],
                        api_key=ncbi_api_key,
                        cache=fetch_cache,
                        doi=paper['doi'],
                        pmcid=paper['pmcid'],
                    )
                    if text:
                        n_pmc += 1

                # ── Unpaywall fallback ─────────────────────────────────────
                if try_uw and paper['doi'] and not text:
                    text = self._fetch_unpaywall_pdf(
                        paper['doi'], unpaywall_email, fetch_cache)
                    if text:
                        n_uw += 1

                if text:
                    paper['full_text']    = text
                    paper['is_full_text'] = True
                    n_pdfs += 1

                if i % 50 == 0 or i == n_abs:
                    print(f"      {i}/{n_abs} processed  "
                          f"(+{n_pmc} PMC, +{n_uw} Unpaywall so far)")

        n_total    = len(self.papers)
        n_abs_only = n_total - n_pdfs
        print(f"    {n_total} papers  |  {n_pdfs} with full-text "
              f"(local PDF + {n_pmc} PMC + {n_uw} Unpaywall)  "
              f"|  {n_abs_only} abstract-only")
        self._loaded = True
        return self

    # ── Corpus search ──────────────────────────────────────────────────────

    def search(self, query_terms: list[str],
               max_sentences_per_paper: int = 12) -> list[dict]:
        """
        Find corpus papers that mention any of the query terms.

        Returns a list of match dicts, each with:
          pmid, title, matched_term, sentences (list[str]), is_full_text
        Full-text papers are returned first.
        """
        if not self._loaded:
            self.load()

        results: list[dict] = []
        qt_lower = [q.lower() for q in query_terms]

        for key, paper in self.papers.items():
            text_src  = paper.get('full_text') or paper.get('abstract', '')
            title_src = paper.get('title', '')
            combined  = f"{title_src}. {text_src}" if text_src else title_src
            if not combined.strip():
                continue

            combined_lower = combined.lower()
            matched_term   = next(
                (q for ql, q in zip(qt_lower, query_terms)
                 if ql in combined_lower),
                None,
            )
            if not matched_term:
                continue

            # Split into sentences, keep those that mention any query term
            raw_sentences = re.split(r'(?<=[.!?])\s+(?=[A-Z])', combined)
            relevant = [
                s.strip()
                for s in raw_sentences
                if any(ql in s.lower() for ql in qt_lower)
                and 25 < len(s.strip()) < 1500
            ][:max_sentences_per_paper]

            if relevant:
                results.append({
                    'pmid':         paper['pmid'],
                    'title':        paper['title'],
                    'matched_term': matched_term,
                    'sentences':    relevant,
                    'is_full_text': paper['is_full_text'],
                })

        # Full-text papers first, then abstract-only
        results.sort(key=lambda r: (0 if r['is_full_text'] else 1))
        return results


# ── Zotero corpus sentiment analysis (Claude) ─────────────────────────────────

def _empty_corpus_result(n_papers: int = 0, n_sentences: int = 0) -> dict:
    return {
        'corpus_papers':          n_papers,
        'corpus_sentences':       n_sentences,
        'corpus_risk_count':      0,
        'corpus_protective_count': 0,
        'corpus_neutral_count':   0,
        'corpus_sentiment_score': 0.0,
        'corpus_summary':         '',
    }


def analyze_corpus_sentiment(
    feature: str,
    phenotype: str,
    paper_matches: list[dict],
    client,
) -> dict:
    """
    Use Claude claude-opus-4-6 (adaptive thinking) to classify the direction of
    association between *feature* and *phenotype* for each relevant sentence
    found in the Zotero corpus.

    Parameters
    ----------
    feature       : query term (gene symbol, rs-ID, HLA allele, clinical measure)
    phenotype     : phenotype string (e.g. "type 2 diabetes")
    paper_matches : output of ZoteroCorpus.search()
    client        : anthropic.Anthropic() instance

    Returns
    -------
    dict with keys: corpus_papers, corpus_sentences, corpus_risk_count,
                    corpus_protective_count, corpus_neutral_count,
                    corpus_sentiment_score (-1 = protective … +1 = risk),
                    corpus_summary
    """
    if not paper_matches or client is None:
        return _empty_corpus_result()

    # Collect sentences, preferring full-text over abstract-only
    all_entries: list[dict] = []
    for match in paper_matches:
        for s in match['sentences']:
            all_entries.append({
                'pmid':        match['pmid'],
                'title':       match['title'][:80],
                'sentence':    s,
                'is_full_text': match['is_full_text'],
            })

    if not all_entries:
        return _empty_corpus_result(len(paper_matches))

    # Cap: prefer full-text sentences up to CORPUS_MAX_SENTENCES
    if len(all_entries) > CORPUS_MAX_SENTENCES:
        ft   = [e for e in all_entries if e['is_full_text']]
        ab   = [e for e in all_entries if not e['is_full_text']]
        n_ft = min(len(ft), int(CORPUS_MAX_SENTENCES * 0.7))
        n_ab = CORPUS_MAX_SENTENCES - n_ft
        all_entries = ft[:n_ft] + ab[:n_ab]

    sentence_block = '\n'.join(
        f"[{i}] (PMID:{e['pmid'] or 'n/a'}, "
        f"{'full-text' if e['is_full_text'] else 'abstract'}) {e['sentence']}"
        for i, e in enumerate(all_entries)
    )

    prompt = (
        f'You are a biomedical expert assessing the association between '
        f'"{feature}" and "{phenotype}".\n\n'
        f'The following {len(all_entries)} sentences come from '
        f'{len(paper_matches)} papers in a curated scientific literature '
        f'corpus. Full-text papers are prioritised over abstracts.\n\n'
        f'For each sentence classify the association as exactly one of:\n'
        f'  RISK        — feature increases risk, severity, or progression\n'
        f'  PROTECTIVE  — feature reduces risk, is beneficial, or protective\n'
        f'  NEUTRAL     — feature mentioned but direction unclear or methodological\n\n'
        f'Return ONLY a JSON object (no markdown, no preamble):\n'
        f'{{\n'
        f'  "classifications": [\n'
        f'    {{"idx": 0, "label": "RISK", "confidence": 0.85}},\n'
        f'    ...\n'
        f'  ],\n'
        f'  "summary": "2-3 sentence synthesis of the overall association"\n'
        f'}}\n\n'
        f'Sentences:\n{sentence_block}'
    )

    try:
        with client.messages.stream(
            model='claude-opus-4-6',
            max_tokens=4096,
            thinking={'type': 'adaptive'},
            messages=[{'role': 'user', 'content': prompt}],
        ) as stream:
            message = stream.get_final_message()

        response_text = ''.join(
            b.text for b in message.content if hasattr(b, 'text')
        )
        # Strip accidental markdown fences
        response_text = re.sub(r'^```[a-z]*\n?', '', response_text.strip(),
                                flags=re.MULTILINE)
        response_text = re.sub(r'\n?```$', '', response_text.strip(),
                                flags=re.MULTILINE)
        result = json.loads(response_text)

        classifications = result.get('classifications', [])
        summary         = result.get('summary', '')

        risk_count  = sum(1 for c in classifications if c.get('label') == 'RISK')
        prot_count  = sum(1 for c in classifications if c.get('label') == 'PROTECTIVE')
        neut_count  = sum(1 for c in classifications if c.get('label') == 'NEUTRAL')
        total       = risk_count + prot_count + neut_count

        # Weighted sentiment: +1 = all RISK, -1 = all PROTECTIVE
        sentiment_score = 0.0
        if total > 0:
            weighted = sum(
                c.get('confidence', 0.5) * (
                    1.0 if c.get('label') == 'RISK'
                    else -1.0 if c.get('label') == 'PROTECTIVE'
                    else 0.0
                )
                for c in classifications
            )
            sentiment_score = round(weighted / total, 4)

        return {
            'corpus_papers':           len(paper_matches),
            'corpus_sentences':        len(all_entries),
            'corpus_risk_count':       risk_count,
            'corpus_protective_count': prot_count,
            'corpus_neutral_count':    neut_count,
            'corpus_sentiment_score':  sentiment_score,
            'corpus_summary':          summary,
        }

    except Exception as exc:
        print(f"    [WARN] Corpus sentiment failed for '{feature}' × "
              f"'{phenotype}': {exc}")
        # corpus_papers=0 ensures compute_combined_score skips the corpus branch
        # so the failure doesn't inject a phantom neutral score or dilute other evidence.
        return _empty_corpus_result(n_papers=0, n_sentences=0)


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
    if any('high_risk' in p for p in parts):
        return 'high_risk'
    if any('low_control' in p for p in parts):
        return 'low_control'
    return 'both'


def _scores_root(pheno_data: str) -> str:
    """
    Return the 'scores/' directory, handling two common layouts:

      Layout A — pheno_data IS the analysis dir:
        {pheno_data}/scores/cohortBnmf/…

      Layout B — pheno_data is the phenotype root (one level above):
        {pheno_data}/combinedAnalysis/scores/cohortBnmf/…

    Prefers Layout A; falls back to Layout B when scores/ is absent.
    """
    direct = os.path.join(pheno_data, 'scores')
    if os.path.isdir(direct):
        return direct
    via_combined = os.path.join(pheno_data, 'combinedAnalysis', 'scores')
    if os.path.isdir(via_combined):
        return via_combined
    # Neither found — return the direct path; caller will produce a clear error
    return direct


def _resolve_cluster_dir(pheno_data, population, analysis_mode, cohort=None):
    """
    Resolve the bNMF output directory from structured path components.

    analysis_mode
    -------------
    cohort        → scores/cohortBnmf/{cohort}/
                    or scores/cohortBnmf/{cohort}/{population}/
                    when the NMF was run with --population separate and a
                    population-specific sub-directory exists.
    combined      → scores/combinedCohortBnmf_{population}/
    clinical_only → scores/combinedCohortBnmf_{population}_clinical_only/
    genomic_only  → scores/combinedCohortBnmf_{population}_genomic_only/

    pheno_data can be either the direct analysis directory (containing scores/)
    or the phenotype root (containing combinedAnalysis/scores/).
    """
    scores = _scores_root(pheno_data)
    if analysis_mode == 'cohort':
        if not cohort:
            raise ValueError("--cohort is required when --analysis_mode cohort")
        cohort_base = os.path.join(scores, 'cohortBnmf', cohort)
        # When the NMF was run with --population separate, outputs land in
        # {cohort}/{population}/ sub-directories.  Check for that first so
        # the correct feature_loadings.csv is found.
        if population not in ('both',):
            pop_subdir = os.path.join(cohort_base, population)
            if os.path.exists(os.path.join(pop_subdir, 'feature_loadings.csv')):
                return pop_subdir
        return cohort_base
    pop_tag = f'_{population}' if population != 'both' else ''
    mode_suffix = {
        'combined':      '',
        'clinical_only': '_clinical_only',
        'genomic_only':  '_genomic_only',
    }.get(analysis_mode, '')
    return os.path.join(scores, f'combinedCohortBnmf{pop_tag}{mode_suffix}')


# Population sub-directory names that should never be treated as cohort names
_POPULATION_DIR_NAMES = {'high_risk', 'low_control', 'both', 'separate'}


def _list_cohort_dirs(pheno_data):
    """
    Return sorted cohort subdirectory names under cohortBnmf/.

    Skips population-label directories (high_risk, low_control) that appear
    at the top level when a cohort was run with population='separate'.
    """
    base = os.path.join(_scores_root(pheno_data), 'cohortBnmf')
    if not os.path.isdir(base):
        return []
    return sorted(
        d for d in os.listdir(base)
        if os.path.isdir(os.path.join(base, d))
        and d not in _POPULATION_DIR_NAMES
    )


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
                           pubtator3_count, clinvar_sig,
                           corpus_papers=0, corpus_sentiment_score=0.0):
    """
    Combine evidence from DisGeNET, Open Targets, PubTator3 and ClinVar
    into a single [0, 1] score.

    PubMed keyword counts are deliberately excluded from the formula here
    (they are stored in the evidence table for reference) because PubTator3
    entity counts already capture publication evidence more precisely.

    When corpus evidence is available (corpus_papers > 0), the Zotero corpus
    sentiment score is incorporated and the base source weights are rescaled
    to occupy 80% of the total, with WEIGHT_CORPUS (20%) allocated to the
    corpus term.

    Weights (no corpus): DisGeNET=0.25, Open Targets=0.30, PubTator3=0.25, ClinVar=0.20
    Weights (corpus):    above × 0.80  +  corpus × 0.20
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

    if corpus_papers > 0:
        # corpus_sentiment_score ∈ [-1, +1]; map to [0, 1] for weighting
        corpus_val = (corpus_sentiment_score + 1.0) / 2.0

        combined = (
            WEIGHT_DISGENET_CORP     * disgenet_score     +
            WEIGHT_OPEN_TARGETS_CORP * open_targets_score +
            WEIGHT_PUBTATOR3_CORP    * pubtator3_norm     +
            WEIGHT_CLINVAR_CORP      * clinvar_val         +
            WEIGHT_CORPUS            * corpus_val
        )
    else:
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


# ── Background data loaders ────────────────────────────────────────────────────

def _find_background_dir(cluster_dir: str, pheno_data: str | None = None) -> str | None:
    """
    Locate the 'background/' directory that holds phenotype_background.json/csv.

    Search order
    ------------
    1. {pheno_data}/background/            — exact match when pheno_data given
    2. {pheno_data}/../background/         — when pheno_data = …/combinedAnalysis/
       and background sits one level up at the phenotype root
    3. Walk upward from cluster_dir (up to 8 levels) looking for a sibling
       'background/' — handles arbitrary nesting without needing pheno_data

    Returns the absolute path to the background directory, or None.
    """
    # ── Direct pheno_data hints (most reliable) ────────────────────────────
    if pheno_data:
        for base in (pheno_data, os.path.dirname(os.path.abspath(pheno_data))):
            candidate = os.path.join(base, 'background')
            if os.path.isdir(candidate):
                return candidate

    # ── Walk up from cluster_dir ───────────────────────────────────────────
    path = os.path.abspath(cluster_dir)
    for _ in range(8):
        parent = os.path.dirname(path)
        candidate = os.path.join(parent, 'background')
        if os.path.isdir(candidate):
            return candidate
        if parent == path:
            break
        path = parent
    return None


def load_background_data(path: str) -> dict:
    """
    Load phenotype background data from a JSON or CSV file produced by
    phenotype_background_analysis.py.

    Returns
    -------
    dict  {term_lower: {'term': str, 'category': str, 'count': int,
                        'sentiment_score': float, 'association_type': str}}

    JSON format  : root['terms'] is a list of dicts with those fields.
    CSV format   : rows with columns category, term, count, sentiment_score,
                   association_type (all optional except term).
    """
    if not path or not os.path.exists(path):
        return {}

    try:
        if path.endswith('.json'):
            import json as _json
            with open(path) as fh:
                raw = _json.load(fh)
            items = raw.get('terms', []) if isinstance(raw, dict) else raw
        else:
            import pandas as _pd
            items = _pd.read_csv(path).to_dict('records')
    except Exception as exc:
        print(f"  [WARN] Could not read background data '{path}': {exc}")
        return {}

    result: dict = {}
    for item in items:
        term = str(item.get('term', '')).strip()
        if not term:
            continue
        result[term.lower()] = {
            'term':             term,
            'category':         str(item.get('category', '')),
            'count':            int(item.get('count', 0) or 0),
            'sentiment_score':  float(item.get('sentiment_score', 0.0) or 0.0),
            'association_type': str(item.get('association_type', '')),
        }

    print(f"  Background data   : {len(result)} terms from {os.path.basename(path)}")
    cats = {}
    for v in result.values():
        cats[v['category']] = cats.get(v['category'], 0) + 1
    for cat, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"    {cat:30s} {n:4d} terms")
    return result


def _background_weight(count: int, sentiment: float) -> float:
    """
    Compute a background evidence weight in [0, 1] that rewards:
      - High publication count   (log-scaled; saturates ~1e6)
      - High sentiment_score     (topic-specificity proxy)

    Formula: log10(count + 2) / 7.0  ×  (0.5 + 0.5 × clamp(sentiment / 0.4))
      - log10 denominator 7 ≈ log10(1e7), so even the most-cited terms < 1
      - sentiment factor ranges from 0.5 (unscored) to 1.0 (highly specific)
    """
    import math
    count_factor     = math.log10(max(count, 0) + 2) / 7.0
    sentiment_clamped = min(1.0, float(sentiment or 0.0) / 0.4)
    sentiment_factor  = 0.5 + 0.5 * sentiment_clamped
    return min(1.0, count_factor * sentiment_factor)


def load_phenotypes_from_background_table(
    background_table_path: str,
    primary_phenotype: str | None = None,
    association_type: str | None = None,
    min_count: int = 0,
) -> list[str]:
    """
    Load all unique terms from a phenotype_background CSV or JSON file produced
    by phenotype_background_analysis.py and return them as a phenotype list for
    use in enrich_with_literature().  Supports both .csv and .json formats.

    Every gene/SNP feature will be queried against every term in this list,
    producing a full (feature × background_term) evidence matrix so that
    each cluster feature has scores for "type 2 diabetes" AND
    "beta-cell dysfunction" AND "metabolic syndrome" etc. simultaneously.

    Parameters
    ----------
    background_table_path : str  — path to phenotype_background.csv or .json
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
    bg_meta = load_background_data(background_table_path)
    if not bg_meta:
        return [primary_phenotype] if primary_phenotype else []

    # Apply filters
    filtered = [
        v for v in bg_meta.values()
        if (not association_type or v['association_type'] == association_type)
        and v['count'] >= min_count
    ]

    # Sort: highest count first (within each category, counts proxy specificity)
    filtered.sort(key=lambda v: -v['count'])

    terms_ordered: list[str] = []
    seen: set[str] = set()
    for v in filtered:
        t = v['term']
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
    auto_background=True,
    pheno_data=None,
    disgenet_api_key=None,
    ncbi_api_key=None,
    unpaywall_email=None,
    top_n=20,
    population_context='auto',
    append_mode=True,
    zotero_corpus_collection=None,
    zotero_db=None,
    anthropic_api_key=None,
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

    When zotero_corpus_collection is provided, full-text and abstract-level
    sentence analysis is performed using a local Zotero library and Claude
    claude-opus-4-6 (adaptive thinking).  The corpus sentiment score is folded into
    the combined evidence score (20% weight; base weights rescaled to 80%).

    Parameters
    ----------
    cluster_dir               : str — bNMF output directory (contains feature_loadings.csv)
    phenotypes                : list[str] — phenotype terms to query; used as the base list
                                 and also as the primary_phenotype anchor when background_table
                                 is supplied (so the primary phenotype is always included)
    background_table          : str | None — path to phenotype_background.csv produced by
                                 phenotype_background_analysis.py.  When supplied, ALL unique
                                 terms in the table are appended to the phenotypes list so
                                 that every feature × every background term is queried.
                                 No ranking or truncation is applied — all terms are used.
    disgenet_api_key          : str | None — DisGeNET API key (free at disgenet.org)
    ncbi_api_key              : str | None — NCBI API key (3→10 req/sec); also used for
                                 PMC full-text fetch of abstract-only Zotero papers
    unpaywall_email           : str | None — email passed to Unpaywall API (required by their
                                 ToS); enables OA PDF fetch for abstract-only Zotero papers
    top_n                     : int — top features per cluster to query
    population_context        : str — 'auto' | 'both' | 'high_risk' | 'low_control'
        Controls which additional queries are run:
          'high_risk'   — standard risk-oriented queries (default behaviour)
          'low_control' — adds protective/prevention queries per feature × phenotype;
                          output columns include protective_count, protective_sentiment,
                          and unique_protective_titles (findings absent from risk queries)
          'both'/'auto' — auto-detect from cluster_dir path; if not detectable, uses 'both'
    append_mode               : bool — if True and output CSVs already exist, new columns are
                                 merged into them rather than overwriting (default: True)
    zotero_corpus_collection  : str | None — Zotero collection name for corpus enrichment
                                 (e.g. 'T2D').  When set, full-text PDFs and abstracts from
                                 the collection are searched and Claude performs sentence-level
                                 sentiment analysis per (feature × phenotype) pair.
    zotero_db                 : str | None — path to zotero.sqlite (default: ~/Zotero/zotero.sqlite)
    anthropic_api_key         : str | None — Anthropic API key for Claude corpus sentiment
                                 analysis; falls back to ANTHROPIC_API_KEY env var if unset

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

    # ── Auto-detect background file if not supplied ────────────────────────
    bg_path = background_table
    if not bg_path and auto_background:
        bg_dir = _find_background_dir(cluster_dir, pheno_data=pheno_data)
        if bg_dir:
            for ext in ('.json', '.csv'):
                candidate = os.path.join(bg_dir, f'phenotype_background{ext}')
                if os.path.exists(candidate):
                    bg_path = candidate
                    print(f"  Background file   : auto-detected {candidate}")
                    break

    # ── Load background metadata (category / count / sentiment) ───────────
    # term_meta maps term_lower → {term, category, count, sentiment_score}
    # Used for: (1) expanding phenotype list; (2) annotating output columns;
    #           (3) computing background_weighted_score.
    term_meta: dict = load_background_data(bg_path) if bg_path else {}

    # ── Expand phenotypes from background data ─────────────────────────────
    if term_meta:
        primary = phenotypes[0] if phenotypes else None
        bg_terms = load_phenotypes_from_background_table(
            bg_path,
            primary_phenotype=primary,
        )
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
    print(f"  Loaded feature loadings: {H_df.shape[0]} clusters × {H_df.shape[1]} features")

    # ── Feature → top-cluster mapping (highest loading across clusters) ────
    # Used to annotate enrichment results with the dominant cluster per feature.
    feat_top_cluster: dict = H_df.idxmax(axis=0).to_dict()

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

    # ── Supplementary features from importantFeaturesAcrossCohortsAndTrainingData ─
    # Load the tiered genomic feature table produced by calculate_top_features_in_cohort.py.
    # Any features NOT already present in feature_loadings.csv are added to the
    # enrichment universe with synthetic loading = normalised matching_z_score and
    # top_cluster = 'supplementary'.  Only HighCases training_data rows are used.
    supp_loadings: dict[str, float] = {}     # feature → synthetic loading (0-1)
    if pheno_data:
        scores_root = _scores_root(pheno_data)
        supp_path = os.path.join(
            scores_root, 'importantCohortScores',
            'importantFeaturesAcrossCohortsAndTrainingData.Filtered.tiered.csv',
        )
        if os.path.exists(supp_path):
            supp_df = pd.read_csv(supp_path)
            # Filter to HighCases training split only
            if 'training_data' in supp_df.columns:
                supp_df = supp_df[supp_df['training_data'] == 'HighCases']
            if 'feature_set' in supp_df.columns and 'cohort' not in supp_df.columns:
                supp_df.rename(columns={'feature_set': 'cohort'}, inplace=True)
            # Normalise feature names: strip gen_ / clin_ prefix for comparison
            nmf_bare = {
                re.sub(r'^(?:gen|clin)_', '', f).lower()
                for f in H_df.columns
            }
            # Build normalised loading from matching_z_score
            z_col = 'matching_z_score' if 'matching_z_score' in supp_df.columns else None
            if z_col and 'feature' in supp_df.columns:
                # Keep best (max z) row per feature across cohorts
                best = (supp_df.groupby('feature')[z_col]
                        .max().reset_index().rename(columns={z_col: 'z'}))
                z_max = best['z'].max() or 1.0
                for _, row in best.iterrows():
                    bare = str(row['feature']).lower()
                    if bare not in nmf_bare:
                        supp_loadings[row['feature']] = round(row['z'] / z_max, 6)
            if supp_loadings:
                print(f"  Supplementary features: {len(supp_loadings)} extra features "
                      f"from importantFeaturesAcrossCohortsAndTrainingData "
                      f"(not in NMF matrix)")
        else:
            print(f"  [INFO] Supplementary features file not found: {supp_path}")

    # Merge supplementary features
    for feat, loading in supp_loadings.items():
        feat_top_cluster[feat] = 'supplementary'
        # Add to all_features; do NOT add to any cluster_top bucket — they are
        # processed separately below using their synthetic loading.
        if feat not in all_features:
            all_features.append(feat)

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

    # ── Load Zotero corpus (optional) ─────────────────────────────────────
    zotero_corpus: ZoteroCorpus | None = None
    anthropic_client = None

    if zotero_corpus_collection:
        if not _ANTHROPIC_AVAILABLE:
            print("  [WARN] --zotero_corpus supplied but 'anthropic' package not "
                  "installed — corpus enrichment disabled.  "
                  "Install with: pip install anthropic")
        else:
            db_path = zotero_db or DEFAULT_ZOTERO_DB
            zotero_corpus = ZoteroCorpus(
                collection_name=zotero_corpus_collection,
                db_path=db_path,
            ).load(ncbi_api_key=ncbi_api_key,
                   unpaywall_email=unpaywall_email)
            if not zotero_corpus.papers:
                print("  [WARN] Zotero corpus is empty — corpus enrichment disabled.")
                zotero_corpus = None
            else:
                api_key = anthropic_api_key or os.environ.get('ANTHROPIC_API_KEY')
                if not api_key:
                    print("  [WARN] No Anthropic API key found — corpus sentiment "
                          "analysis disabled.  Set ANTHROPIC_API_KEY or pass "
                          "--anthropic_api_key.")
                    zotero_corpus = None
                else:
                    anthropic_client = _anthropic.Anthropic(api_key=api_key)
                    print(f"  Corpus sentiment: enabled "
                          f"({len(zotero_corpus.papers)} papers via Claude claude-opus-4-6)")

    corpus_active = zotero_corpus is not None and anthropic_client is not None

    sources_label = "DisGeNET, Open Targets, PubTator3, PubMed, ClinVar"
    if corpus_active:
        sources_label += f", Zotero corpus ({zotero_corpus_collection})"
    print(f"\n  Evidence sources: {sources_label}")
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

        # ── Zotero corpus sentiment (optional) ────────────────────────────
        corpus_result = _empty_corpus_result()
        if corpus_active:
            paper_matches = zotero_corpus.search([query_term, term])
            if paper_matches:
                corpus_result = analyze_corpus_sentiment(
                    feature=query_term,
                    phenotype=pheno,
                    paper_matches=paper_matches,
                    client=anthropic_client,
                )

        combined = compute_combined_score(
            disgenet_score, ot_score, pt3_count, clinvar_sig,
            corpus_papers=corpus_result['corpus_papers'],
            corpus_sentiment_score=corpus_result['corpus_sentiment_score'],
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

        # ── Background category / weight annotation ────────────────────────
        pheno_meta  = term_meta.get(pheno.lower(), {})
        bg_category = pheno_meta.get('category', '')
        bg_count    = pheno_meta.get('count', 0)
        bg_sentiment = pheno_meta.get('sentiment_score', 0.0)
        bg_weight   = _background_weight(bg_count, bg_sentiment)

        row = {
            'term':                  term,
            'query_term':            query_term,
            'phenotype':             pheno,
            'background_category':   bg_category,
            'background_count':      bg_count,
            'background_sentiment':  round(bg_sentiment, 4),
            'background_weight':     round(bg_weight, 4),
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
        }
        if corpus_active:
            row.update(corpus_result)
        evidence_rows.append(row)

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
                    score      = ev_lookup.get((term, pheno), 0.0)
                    pheno_meta = term_meta.get(pheno.lower(), {})
                    bg_cat     = pheno_meta.get('category', '')
                    bg_wt      = _background_weight(
                        pheno_meta.get('count', 0),
                        pheno_meta.get('sentiment_score', 0.0),
                    )
                    cluster_ranked_rows.append({
                        'cluster':                   cluster,
                        'feature':                   feat,
                        'top_cluster':               feat_top_cluster.get(feat, cluster),
                        'term':                      term,
                        'phenotype':                 pheno,
                        'background_category':       bg_cat,
                        'background_weight':         round(bg_wt, 4),
                        'cluster_rank':              rank,
                        'loading':                   round(loading, 6),
                        'combined_score':            score,
                        'weighted_rank_score':       round(loading * score, 6),
                        'background_weighted_score': round(loading * score * bg_wt, 6),
                    })

    # ── Supplementary features (not in NMF; from importantFeaturesAcrossCohortsAndTrainingData) ─
    for rank, (feat, loading) in enumerate(
            sorted(supp_loadings.items(), key=lambda x: x[1], reverse=True), start=1):
        for term in all_terms.get(feat, []):
            for pheno in phenotypes:
                score      = ev_lookup.get((term, pheno), 0.0)
                pheno_meta = term_meta.get(pheno.lower(), {})
                bg_cat     = pheno_meta.get('category', '')
                bg_wt      = _background_weight(
                    pheno_meta.get('count', 0),
                    pheno_meta.get('sentiment_score', 0.0),
                )
                cluster_ranked_rows.append({
                    'cluster':                   'supplementary',
                    'feature':                   feat,
                    'top_cluster':               'supplementary',
                    'term':                      term,
                    'phenotype':                 pheno,
                    'background_category':       bg_cat,
                    'background_weight':         round(bg_wt, 4),
                    'cluster_rank':              rank,
                    'loading':                   loading,
                    'combined_score':            score,
                    'weighted_rank_score':       round(loading * score, 6),
                    'background_weighted_score': round(loading * score * bg_wt, 6),
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
    # Use background_weighted_score when available, else weighted_rank_score
    summary_sort_col = (
        'background_weighted_score'
        if 'background_weighted_score' in cluster_ranked_df.columns
        and cluster_ranked_df['background_weighted_score'].sum() > 0
        else 'weighted_rank_score'
    )
    summary_rows = []
    if not cluster_ranked_df.empty:
        for cluster, grp in cluster_ranked_df.groupby('cluster'):
            best = grp.nlargest(1, summary_sort_col).iloc[0]
            # Top hit per background category (when categories available)
            category_bests = {}
            if 'background_category' in grp.columns:
                for cat, cat_grp in grp[grp['background_category'] != ''].groupby(
                        'background_category'):
                    cb = cat_grp.nlargest(1, summary_sort_col).iloc[0]
                    category_bests[cat] = (cb['term'], round(cb[summary_sort_col], 6))

            row = {
                'cluster':             cluster,
                'top_term':            best['term'],
                'top_feature':         best['feature'],
                'top_cluster':         best.get('top_cluster', cluster),
                'dominant_phenotype':  best['phenotype'],
                'dominant_category':   best.get('background_category', ''),
                'max_evidence_score':  best['combined_score'],
                'weighted_rank_score': best['weighted_rank_score'],
            }
            if 'background_weighted_score' in best.index:
                row['background_weighted_score'] = best['background_weighted_score']
            # Add per-category top terms as flat columns
            for cat, (cat_term, cat_score) in sorted(category_bests.items()):
                row[f'top_{cat}_term']  = cat_term
                row[f'top_{cat}_score'] = cat_score
            summary_rows.append(row)

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
            cat_tag = (f"  [{row['dominant_category']}]"
                       if row.get('dominant_category') else '')
            bw_tag  = (f"  bg_wt={row['background_weighted_score']:.4f}"
                       if 'background_weighted_score' in row.index else '')
            print(f"    {row['cluster']:20s} → {row['top_term']} "
                  f"× {row['dominant_phenotype']}{cat_tag}  "
                  f"(score={row['max_evidence_score']:.3f}{bw_tag})")

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
        '--cluster_dir', default=None,
        help=(
            "Path to bNMF output directory containing feature_loadings.csv.\n"
            "If omitted, the path is auto-resolved from --pheno_data,\n"
            "--population, --analysis_mode (and --cohort for cohort mode)."
        ),
    )
    parser.add_argument(
        '--pheno_data', default=None,
        help=(
            "Base directory of the analysis (e.g. results/type2Diabetes/combinedAnalysis).\n"
            "Used together with --population and --analysis_mode to auto-resolve\n"
            "--cluster_dir when it is not provided explicitly."
        ),
    )
    parser.add_argument(
        '--population',
        default='high_risk',
        choices=['high_risk', 'low_control', 'both'],
        help=(
            "Population group to analyse (default: high_risk).\n"
            "Controls which combinedCohortBnmf sub-directory is used when\n"
            "--cluster_dir is auto-resolved."
        ),
    )
    parser.add_argument(
        '--analysis_mode',
        default='cohort',
        choices=['cohort', 'combined', 'clinical_only', 'genomic_only'],
        help=(
            "Which bNMF results to use (default: cohort):\n"
            "  cohort        — per-cohort results under cohortBnmf/{cohort}/\n"
            "  combined      — cross-cohort results (combinedCohortBnmf_{population})\n"
            "  clinical_only — clinical-only NMF (combinedCohortBnmf_{population}_clinical_only)\n"
            "  genomic_only  — genomic-only NMF (combinedCohortBnmf_{population}_genomic_only)"
        ),
    )
    parser.add_argument(
        '--cohort', default=None,
        help=(
            "Cohort name for --analysis_mode cohort (e.g. 'cardio').\n"
            "If omitted, all subdirectories under cohortBnmf/ are processed."
        ),
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
            "Path to phenotype_background.csv or .json produced by\n"
            "phenotype_background_analysis.py.  When supplied, all unique\n"
            "terms in the table (subtypes, co-phenotypes, mechanisms, pathways,\n"
            "complications, comorbidities) are added to the phenotype list so\n"
            "every gene/SNP feature is scored against every background term.\n"
            "If omitted, a phenotype_background.json/.csv is auto-detected from\n"
            "a 'background/' directory adjacent to the analysis root.\n"
            "Outputs gain: background_category, background_count,\n"
            "background_sentiment, background_weight, background_weighted_score."
        ),
    )
    parser.add_argument(
        '--no_background_auto', action='store_true',
        help=(
            "Disable auto-detection of the phenotype_background file.\n"
            "Use when you want to run without background terms even if a\n"
            "background/ directory exists alongside the analysis."
        ),
    )
    parser.add_argument(
        '--top_n', type=int, default=20,
        help="Top features per cluster to query (default: 20)",
    )
    parser.add_argument(
        '--ncbi_api_key', default=None,
        help=(
            "NCBI Entrez API key — increases rate limit from 3 to 10 req/sec.\n"
            "Also used for PubMed Central full-text fetch of abstract-only "
            "Zotero papers."
        ),
    )
    parser.add_argument(
        '--unpaywall_email', default=None,
        metavar='EMAIL',
        help=(
            "Email address to pass to the Unpaywall API (required by their ToS).\n"
            "When set, abstract-only Zotero papers are checked against Unpaywall\n"
            "for legal open-access PDFs; those found are downloaded and extracted\n"
            "as full text.  Runs after PMC fetch (PMC takes priority)."
        ),
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
    parser.add_argument(
        '--zotero_corpus', default=None,
        metavar='COLLECTION',
        help=(
            "Zotero collection name to use as a curated literature corpus for\n"
            "sentence-level sentiment analysis (e.g. 'T2D').  Full-text PDFs are\n"
            "prioritised over abstracts.  Requires the 'anthropic' Python package\n"
            "and a valid Anthropic API key (--anthropic_api_key or ANTHROPIC_API_KEY)."
        ),
    )
    parser.add_argument(
        '--zotero_db', default=None,
        metavar='PATH',
        help=(
            f"Path to Zotero SQLite database (default: {DEFAULT_ZOTERO_DB}).\n"
            "Only used when --zotero_corpus is set."
        ),
    )
    parser.add_argument(
        '--anthropic_api_key', default=None,
        help=(
            "Anthropic API key for Claude corpus sentiment analysis.\n"
            "Falls back to the ANTHROPIC_API_KEY environment variable if not set.\n"
            "Only used when --zotero_corpus is set."
        ),
    )

    args = parser.parse_args()

    # ── Resolve cluster_dir(s) ────────────────────────────────────────────────
    if args.cluster_dir:
        cluster_dirs = [args.cluster_dir]
    elif args.pheno_data:
        if args.analysis_mode == 'cohort':
            if args.cohort:
                cluster_dirs = [
                    _resolve_cluster_dir(args.pheno_data, args.population,
                                         'cohort', args.cohort)
                ]
            else:
                # No cohort specified → run all cohorts found on disk
                cohorts = _list_cohort_dirs(args.pheno_data)
                if not cohorts:
                    print(f"[ERROR] No cohort sub-directories found under "
                          f"{args.pheno_data}/scores/cohortBnmf/")
                    raise SystemExit(1)
                print(f"  Found {len(cohorts)} cohort(s): {', '.join(cohorts)}")
                cluster_dirs = [
                    _resolve_cluster_dir(args.pheno_data, args.population,
                                         'cohort', c)
                    for c in cohorts
                ]
        else:
            cluster_dirs = [
                _resolve_cluster_dir(args.pheno_data, args.population,
                                      args.analysis_mode)
            ]
    else:
        print("[ERROR] Provide either --cluster_dir or --pheno_data.")
        raise SystemExit(1)

    phenotypes = [p.strip() for p in args.phenotypes.split(',')]
    shared_kw  = dict(
        phenotypes=phenotypes,
        background_table=args.background_table,
        auto_background=not args.no_background_auto,
        pheno_data=args.pheno_data,          # passed to _find_background_dir + supplementary features
        disgenet_api_key=args.disgenet_api_key,
        ncbi_api_key=args.ncbi_api_key,
        unpaywall_email=args.unpaywall_email,
        top_n=args.top_n,
        population_context=args.population_context,
        append_mode=not args.no_append,
        zotero_corpus_collection=args.zotero_corpus,
        zotero_db=args.zotero_db,
        anthropic_api_key=args.anthropic_api_key,
    )

    for cdir in cluster_dirs:
        if not os.path.isdir(cdir):
            print(f"  [SKIP] Directory not found: {cdir}")
            continue
        print(f"\n{'='*70}\n  cluster_dir: {cdir}\n{'='*70}")
        enrich_with_literature(cluster_dir=cdir, **shared_kw)
