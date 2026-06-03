"""
snpToGene.py
------------
Two-stage pipeline driven by a penalized-model filter file.

  STAGE 1 – Build feature→gene lookup
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Features are sourced from filter_file
  (featureScoresReducedFinalModel.combined.filtered.csv).
  Only rows with non-zero 'coefs' are queried.

  Three classes of queryable features are extracted:

    rsID features   (G / GxG)   : 'rs12345' or 'rs12345,rs67890'
    HLA allele features          : 'DRB_5301' → queried as 'HLA-DRB*5301'
    Covariate features           : rows where model == 'covariate' — skipped
    Clinical / other features    : len(feature.split(',')) > 2, or no rsID
                                   and not HLA pattern — skipped

  Previously queried SNPs are loaded from the output directory and skipped,
  allowing incremental re-runs without redundant Entrez calls.

  Queries NCBI Entrez via three tiers (rsID) or gene-db direct search (HLA):

      g_features_snp_map.csv   – one row per rsID with full Entrez summary
      gene_records.csv         – one row per (feature × gene) pair
      unmatched_snps.csv       – rsIDs / HLA features with no gene hit

  STAGE 2 – Annotate feature scores
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Merge gene annotations back onto every row of filter_file.
  Non-queryable rows (clinical, covariate) receive NaN in gene columns.

  Output:
      feature_scores_annotated.csv – original feature_scores with added columns:
        snp1                : rsID or 'HLA-DRB*5301' form
        snp1_gene_names     : semicolon-joined gene names for snp1 (best tier first)
        snp1_gene_ids       : semicolon-joined NCBI gene IDs for snp1
        snp1_gene_source    : 'snp_db' | 'elink' | 'gene_db' | 'hla_gene_db'
        snp2                : second rsID (GxG only)
        snp2_gene_names / snp2_gene_ids / snp2_gene_source  (GxG only)
        entrez_matched      : True if at least one gene was found for this row

Usage
-----
python snpToGene.py \\
    --filter_file  path/to/featureScoresReducedFinalModel.combined.filtered.csv \\
    --email        you@institution.edu \\
    --output_path  ./snpToGene_output
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import re
import logging
import argparse
from pathlib import Path

import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from Bio import Entrez

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DB           = "snp"
DB_GENE      = "gene"
PARAM_EUTILS = {"usehistory": "Y", "retmode": "xml", "retmax": 10}
# Gene db calls do not use web history; XML retmode is still appropriate.
PARAM_GENE   = {"retmode": "xml", "retmax": 10}

_RS_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)

# Gene source tiers — lower number = more authoritative
GENE_SOURCE_RANK = {"snp_db": 1, "elink": 2, "gene_db": 3, "hla_gene_db": 4}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _is_valid_rsid(token: str) -> bool:
    """True if *token* is a bare rsID (rs + digits only)."""
    return bool(_RS_PATTERN.match(token.strip()))


def _rsids_from_feature(feature: str) -> list[str]:
    """
    Return the ordered list of rsIDs encoded in one 'feature' cell.

    Handles:
      'rs1234567'          → ['rs1234567']
      'rs1234567,rs891011' → ['rs1234567', 'rs891011']
      'age' / 'HLA-A' / … → []
    """
    tokens = [t.strip() for t in str(feature).split(",")]
    return [t.lower() for t in tokens if _is_valid_rsid(t)]


def _feature_type(rsids: list[str]) -> str:
    """Classify a feature row based on how many rsIDs it contains."""
    n = len(rsids)
    if n == 1:
        return "G"
    if n >= 2:
        return "GxG"
    return "other"


def _is_hla_feature(feature: str) -> bool:
    """
    Return True if *feature* looks like an HLA allele encoded with '_'.

    Criteria
    --------
    - Not a valid rsID
    - No commas (not a GxG or multi-value clinical feature)
    - Exactly one underscore separator
    - Both parts are non-empty and alphanumeric

    Examples
    --------
    'DRB_5301'  → True   ('HLA-DRB*5301')
    'DRB5_301'  → True   ('HLA-DRB5*301')
    'rs1234'    → False  (rsID)
    'bmi'       → False  (no underscore)
    'blood_pressure_medication' → False  (more than one underscore)
    """
    feature = str(feature).strip()
    if _is_valid_rsid(feature) or ',' in feature:
        return False
    parts = feature.split('_')
    return len(parts) == 2 and all(p and p.isalnum() for p in parts)


def _hla_feature_to_query(feature: str) -> str:
    """
    Convert HLA allele underscore notation to the standard query form.

    'DRB_5301'  → 'HLA-DRB*5301'
    'DRB5_301'  → 'HLA-DRB5*301'

    The formula is: 'HLA-' + split[0] + '*' + split[1]
    """
    parts = str(feature).strip().split('_', 1)
    return f"HLA-{parts[0]}*{parts[1]}"


# ---------------------------------------------------------------------------
# STAGE 1A – Extract queryable features from filter_file
# ---------------------------------------------------------------------------

def extract_features_from_filter(
    filter_file: str,
) -> tuple[set[str], dict[str, str]]:
    """
    Extract queryable rsIDs and HLA features from filter_file.

    Only rows with non-zero 'coefs' are considered.
    Covariate rows (model contains 'covariate') and intercept rows
    (feature == '(Intercept)') are always excluded.

    Classification rules (applied in order)
    ----------------------------------------
    1. Skip features where model contains 'covariate' (case-insensitive).
    2. Skip features where feature == '(Intercept)'.
    3. Skip features where ``coefs == 0``.
    4. Skip features where ``len(feature.split(',')) > 2``.
    5. Try to parse rsIDs via ``_rsids_from_feature``.  Valid hits → snps.
    6. If no rsIDs found and the feature matches the HLA allele pattern → hla_map.
    7. All remaining features are clinical / non-genomic — skipped.

    Returns
    -------
    snps    : set[str]         lower-cased rsIDs to query via SNP db tiers
    hla_map : dict[str, str]   {original_feature_lower: 'HLA-DRB*5301'}
    """
    log.info("Loading filter file: %s", filter_file)
    df = pd.read_csv(filter_file)
    log.info("  %d rows × %d columns loaded.", *df.shape)

    if "feature" not in df.columns:
        raise ValueError(
            f"filter_file must contain a 'feature' column. "
            f"Found: {df.columns.tolist()}"
        )
    if "coefs" not in df.columns:
        raise ValueError(
            f"filter_file must contain a 'coefs' column. "
            f"Found: {df.columns.tolist()}"
        )

    # Identify covariate feature names from any covariate model row, then
    # exclude those feature names across all models.
    covariate_feature_names: set[str] = set()
    if "model" in df.columns:
        cov_model_mask = df["model"].str.contains("covariate", case=False, na=False)
        covariate_feature_names = set(
            df.loc[cov_model_mask, "feature"].dropna().str.strip().unique()
        )
        log.info(
            "  Identified %d covariate feature names (from %d covariate model rows).",
            len(covariate_feature_names), int(cov_model_mask.sum()),
        )

    # Exclude: covariate features in any model + explicit intercept rows
    exclude_mask = (
        df["feature"].str.strip().isin(covariate_feature_names)
        | (df["feature"].str.strip() == "(Intercept)")
    )
    log.info("  Excluding %d rows (covariate features or intercept, across all models).", int(exclude_mask.sum()))

    nonzero_mask = (df["coefs"] != 0) & ~exclude_mask
    log.info(
        "  Keeping %d non-zero coef rows (dropping %d zero-coef rows).",
        int(nonzero_mask.sum()), int((df["coefs"] == 0).sum()),
    )
    candidates = df.loc[nonzero_mask, "feature"].dropna().str.strip().unique()

    snps:    set[str]       = set()
    hla_map: dict[str, str] = {}
    n_clinical = 0

    for feat in candidates:
        feat_lower = feat.lower()

        if len(feat.split(",")) > 2:
            n_clinical += 1
            continue

        rsids = _rsids_from_feature(feat)
        if rsids:
            snps.update(rsids)
            continue

        if _is_hla_feature(feat_lower):
            hla_map[feat_lower] = _hla_feature_to_query(feat_lower)
            continue

        n_clinical += 1

    log.info(
        "  Extracted %d unique rsIDs, %d HLA features; %d clinical/other skipped.",
        len(snps), len(hla_map), n_clinical,
    )
    return snps, hla_map


# ---------------------------------------------------------------------------
# STAGE 1A (cache) – load previously queried results
# ---------------------------------------------------------------------------

def load_cache(
    output_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """
    Load existing Entrez results written by a prior run.

    Reads gene_records.csv, g_features_snp_map.csv, and unmatched_snps.csv
    from *output_path*.  Missing files are silently ignored.

    Returns
    -------
    cached_gene_records : pd.DataFrame   (may be empty)
    cached_snp_records  : pd.DataFrame   (may be empty)
    already_queried     : set[str]       lower-cased rsIDs + HLA keys already done
    """
    already_queried: set[str] = set()

    def _add_to_queried(values: "pd.Series") -> None:
        """Extract individual rsIDs from cache values; keep non-rsID tokens as-is."""
        for val in values.dropna().str.lower().unique():
            rsids = _rsids_from_feature(val)
            if rsids:
                already_queried.update(rsids)
            else:
                already_queried.add(val)   # HLA feature keys, etc.

    # gene_records: snpId column holds the queried rsID / HLA feature key
    gene_p = output_path / "gene_records.csv"
    cached_gene_records = pd.DataFrame()
    if gene_p.exists():
        cached_gene_records = pd.read_csv(gene_p)
        if "snpId" in cached_gene_records.columns:
            _add_to_queried(cached_gene_records["snpId"])
        log.info(
            "Cache: loaded %d gene records from %s",
            len(cached_gene_records), gene_p,
        )

    # unmatched_snps: rsIDs that were queried but returned no gene
    unmatched_p = output_path / "unmatched_snps.csv"
    if unmatched_p.exists():
        df_unmatched = pd.read_csv(unmatched_p)
        if "rsid" in df_unmatched.columns:
            _add_to_queried(df_unmatched["rsid"])
        log.info(
            "Cache: loaded %d unmatched rsIDs from %s",
            len(df_unmatched), unmatched_p,
        )

    # g_features_snp_map: full SNP db summaries
    snp_p = output_path / "g_features_snp_map.csv"
    cached_snp_records = pd.DataFrame()
    if snp_p.exists():
        cached_snp_records = pd.read_csv(snp_p)
        log.info(
            "Cache: loaded %d SNP summary records from %s",
            len(cached_snp_records), snp_p,
        )

    log.info("Cache: %d rsIDs / HLA keys already queried — will skip these.",
             len(already_queried))
    return cached_gene_records, cached_snp_records, already_queried


# ---------------------------------------------------------------------------
# STAGE 1B – Query Entrez  (three-tier gene lookup)
# ---------------------------------------------------------------------------
#
# Tier 1 – SNP db GENES field  (most authoritative; gene_source = 'snp_db')
# Tier 2 – elink SNP→gene      (gene_source = 'elink')
# Tier 3 – gene db direct search (gene_source = 'gene_db')
# ---------------------------------------------------------------------------

def _get_snp_uid(snp_id: str) -> list[str]:
    """Return the NCBI UID list for *snp_id* from the SNP database."""
    handle = Entrez.esearch(db=DB, term=snp_id, **PARAM_EUTILS)
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])


def _parse_esummary(dict_obj: dict) -> dict:
    """Extract the first DocumentSummary from a SNP esummary result; return {} on failure."""
    try:
        return dict_obj["DocumentSummarySet"]["DocumentSummary"][0]
    except (KeyError, IndexError) as exc:
        log.warning("Could not parse SNP eSummary: %s", exc)
        return {}


def _elink_snp_to_gene_ids(snp_uid: str) -> list[str]:
    """
    Return gene db UIDs linked to *snp_uid* via NCBI elink (SNP→gene).

    Uses the 'snp_gene' link name; falls back to any LinkSetDb entry whose
    DbTo is 'gene' in case the link name varies.
    """
    handle = Entrez.elink(dbfrom=DB, db=DB_GENE, id=snp_uid)
    result = Entrez.read(handle)
    handle.close()
    gene_ids: list[str] = []
    for linkset in result:
        for lsdb in linkset.get("LinkSetDb", []):
            if lsdb.get("DbTo") == DB_GENE:
                gene_ids.extend(link["Id"] for link in lsdb.get("Link", []))
    return gene_ids


def _fetch_gene_summaries(gene_ids: list[str], snp_id: str,
                           source: str) -> list[dict]:
    """
    Fetch gene db eSummary records for *gene_ids* and format as gene rows.

    Parameters
    ----------
    gene_ids : list of gene db UIDs (strings)
    snp_id   : rsID these genes are associated with
    source   : 'elink' or 'gene_db'

    Returns
    -------
    list of dicts with keys: GENE_ID, NAME, snpId, gene_source, gene_rank
    """
    if not gene_ids:
        return []
    handle = Entrez.esummary(db=DB_GENE, id=",".join(gene_ids), **PARAM_GENE)
    result = Entrez.read(handle)
    handle.close()

    rows: list[dict] = []
    docs = result.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    for rank, doc in enumerate(docs, start=1):
        uid  = doc.attributes.get("uid", "")
        name = str(doc.get("Name", "")).strip()
        if not name:
            name = str(doc.get("NomenclatureSymbol", "")).strip()
        rows.append({
            "GENE_ID":     uid,
            "NAME":        name,
            "snpId":       snp_id,
            "gene_source": source,
            "gene_rank":   rank,
        })
    return rows


def _search_gene_db(queries: list[str], feature_id: str,
                    source: str, max_hits: int = 5) -> list[dict]:
    """
    Try each query string against the gene database in order; return results
    from the first query that yields hits.
    """
    for term in queries:
        try:
            handle  = Entrez.esearch(db=DB_GENE, term=term, **PARAM_GENE)
            result  = Entrez.read(handle)
            handle.close()
            id_list = result.get("IdList", [])
        except Exception as exc:
            log.warning("  gene db esearch (%r) failed: %s", term, exc)
            id_list = []
        if id_list:
            log.debug("  gene db search %r → %d hit(s)", term, len(id_list))
            return _fetch_gene_summaries(id_list[:max_hits], feature_id, source=source)
    return []


def _search_gene_db_for_snp(snp_id: str) -> list[dict]:
    """Tier 3: search gene db for an rsID not found in the SNP db."""
    return _search_gene_db(
        queries=[f"{snp_id}[SNP]", snp_id],
        feature_id=snp_id,
        source="gene_db",
    )


def _search_gene_db_for_hla(orig_feat: str, hla_query: str) -> list[dict]:
    """
    Search the gene database for an HLA allele feature.

    Tries three progressively broader queries, returning the first batch that
    contains at least one HLA-region gene (NAME starts with 'HLA').

    Query order:
      1. ``'HLA-DRB*5301[gene name]'`` – field-tagged gene name search
      2. ``'HLA-DRB*5301'``            – free-text
      3. ``'HLA-DRB'``                 – gene symbol prefix (fallback)
    """
    gene_prefix = hla_query.split('*')[0]   # e.g. 'HLA-DRB'
    for term in [f"{hla_query}[gene name]", hla_query, gene_prefix]:
        rows = _search_gene_db([term], feature_id=orig_feat, source="hla_gene_db")
        hla_rows = [r for r in rows if str(r.get("NAME", "")).upper().startswith("HLA")]
        if hla_rows:
            return hla_rows
    return []


def query_hla_genes(
    hla_map: dict[str, str],
    skip_keys: set[str] | None = None,
) -> pd.DataFrame:
    """
    Query the NCBI gene database for each HLA allele feature.

    Parameters
    ----------
    hla_map   : {original_feature_lower: 'HLA-DRB*5301', ...}
    skip_keys : lower-cased feature keys already in the cache — skipped

    Returns
    -------
    pd.DataFrame with columns: GENE_ID, NAME, snpId, gene_source, gene_rank
    """
    if not hla_map:
        return pd.DataFrame()

    skip = {k.lower() for k in (skip_keys or set())}
    to_query = {k: v for k, v in hla_map.items() if k not in skip}
    skipped  = len(hla_map) - len(to_query)
    if skipped:
        log.info("Cache: skipping %d already-queried HLA feature(s).", skipped)

    rows: list[dict] = []
    total = len(to_query)
    for i, (orig_feat, hla_query) in enumerate(sorted(to_query.items()), 1):
        log.info("[%d/%d] HLA query: %s → %s", i, total, orig_feat, hla_query)
        try:
            gene_rows = _search_gene_db_for_hla(orig_feat, hla_query)
            if gene_rows:
                rows.extend(gene_rows)
                log.debug("  Found %d gene(s).", len(gene_rows))
            else:
                log.debug("  No genes found for %s.", hla_query)
        except Exception as exc:
            log.warning("  HLA gene query failed for %s: %s", hla_query, exc)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def query_entrez_for_genes(
    snps: set[str],
    skip_snps: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """
    Query NCBI Entrez for gene associations for each rsID using a three-tier
    strategy that maximises recall while preserving source provenance.

    Tier 1 – SNP db GENES field  (gene_source = 'snp_db')
    Tier 2 – elink SNP→gene      (gene_source = 'elink')
    Tier 3 – gene db search      (gene_source = 'gene_db')

    Parameters
    ----------
    snps      : set of lower-cased rsIDs to query
    skip_snps : rsIDs already in the cache — skipped entirely

    Returns
    -------
    snp_records : pd.DataFrame   Full SNP db DocumentSummary per rsID found
    gene_records : pd.DataFrame  One row per (rsID × gene)
    unmatched : set[str]         rsIDs for which no gene was found via any tier
    """
    skip = {s.lower() for s in (skip_snps or set())}
    to_query = sorted(s for s in snps if s not in skip)
    skipped  = len(snps) - len(to_query)
    if skipped:
        log.info("Cache: skipping %d already-queried SNP(s).", skipped)

    snp_rows:  list[dict] = []
    gene_rows: list[dict] = []
    unmatched: set[str]   = set()

    total = len(to_query)
    for i, snp in enumerate(to_query, 1):
        log.info("[%d/%d] %s …", i, total, snp)
        found_genes = False

        # ------------------------------------------------------------------
        # Tier 1: SNP database
        # ------------------------------------------------------------------
        snp_uid = None
        record  = {}
        try:
            id_list = _get_snp_uid(snp)
            snp_uid = id_list[0] if id_list else None
        except Exception as exc:
            log.warning("  esearch(snp) failed: %s", exc)
            id_list = []

        if snp_uid:
            try:
                handle       = Entrez.esummary(db=DB, id=snp_uid, **PARAM_EUTILS)
                summary_dict = Entrez.read(handle)
                handle.close()
                record       = _parse_esummary(summary_dict)
            except Exception as exc:
                log.warning("  esummary(snp, uid=%s) failed: %s", snp_uid, exc)

        if record:
            row           = dict(record)
            row["snpId"]  = snp
            row["IDList"] = id_list
            snp_rows.append(row)

            genes = record.get("GENES", [])
            if genes:
                for rank, g in enumerate(genes, start=1):
                    gr               = dict(g)
                    gr["snpId"]      = snp
                    gr["gene_source"] = "snp_db"
                    gr["gene_rank"]  = rank
                    gene_rows.append(gr)
                log.debug("  Tier 1: %d gene(s) from SNP db GENES field.", len(genes))
                found_genes = True

        # ------------------------------------------------------------------
        # Tier 2: elink SNP→gene  (SNP in db but GENES empty)
        # ------------------------------------------------------------------
        if not found_genes and snp_uid:
            try:
                elink_ids = _elink_snp_to_gene_ids(snp_uid)
            except Exception as exc:
                log.warning("  elink(snp→gene) failed: %s", exc)
                elink_ids = []

            if elink_ids:
                try:
                    rows = _fetch_gene_summaries(elink_ids[:10], snp, source="elink")
                    if rows:
                        gene_rows.extend(rows)
                        log.debug("  Tier 2: %d gene(s) via elink.", len(rows))
                        found_genes = True
                except Exception as exc:
                    log.warning("  elink gene fetch failed: %s", exc)

        # ------------------------------------------------------------------
        # Tier 3: gene database direct search
        # ------------------------------------------------------------------
        if not found_genes:
            try:
                rows = _search_gene_db_for_snp(snp)
                if rows:
                    gene_rows.extend(rows)
                    log.debug("  Tier 3: %d gene(s) from gene db search.", len(rows))
                    found_genes = True
            except Exception as exc:
                log.warning("  gene db search failed: %s", exc)

        if not found_genes:
            log.debug("  No genes found via any tier.")
            unmatched.add(snp)

    snp_records  = pd.DataFrame(snp_rows)
    gene_records = pd.DataFrame(gene_rows)

    tier_counts = (
        gene_records["gene_source"].value_counts().to_dict()
        if not gene_records.empty and "gene_source" in gene_records.columns
        else {}
    )
    log.info(
        "Entrez complete → %d SNP records | %d gene rows "
        "(snp_db=%d, elink=%d, gene_db=%d) | %d unmatched",
        len(snp_rows), len(gene_rows),
        tier_counts.get("snp_db", 0),
        tier_counts.get("elink",  0),
        tier_counts.get("gene_db", 0),
        len(unmatched),
    )
    return snp_records, gene_records, unmatched


# ---------------------------------------------------------------------------
# STAGE 1C – Build the rsID → gene lookup dict
# ---------------------------------------------------------------------------

def build_snp_gene_map(gene_records: pd.DataFrame) -> dict[str, dict]:
    """
    Collapse the per-(rsID × gene) DataFrame into a dict keyed by lower-case rsID.

    When a SNP maps to multiple genes, names and IDs are joined with ';' in
    gene_rank order (most relevant gene first within each tier).

    The most authoritative tier for each SNP is recorded as ``gene_source``
    in the returned dict (snp_db > elink > gene_db).

    Returns
    -------
    {
      'rs1234': {'gene_names': 'BRCA1',     'gene_ids': '672',
                 'gene_source': 'snp_db'},
      'rs5678': {'gene_names': 'TP53;MDM2', 'gene_ids': '7157;4193',
                 'gene_source': 'elink'},
      ...
    }
    """
    if gene_records.empty:
        log.warning("gene_records is empty — snp_gene_map will be empty.")
        return {}

    # Normalise column names (Entrez returns mixed case; upper for safety)
    gr = gene_records.copy()
    gr.columns = [c.upper() for c in gr.columns]

    # After upper-casing: snpId→SNPID, gene_source→GENE_SOURCE, gene_rank→GENE_RANK
    name_col   = "NAME"        if "NAME"        in gr.columns else None
    id_col     = "GENE_ID"     if "GENE_ID"     in gr.columns else None
    source_col = "GENE_SOURCE" if "GENE_SOURCE" in gr.columns else None
    rank_col   = "GENE_RANK"   if "GENE_RANK"   in gr.columns else None

    result: dict[str, dict] = {}
    for snp_id, grp in gr.groupby("SNPID"):
        if source_col and rank_col:
            grp = grp.copy()
            grp["_tier"] = grp[source_col].map(GENE_SOURCE_RANK).fillna(99)
            grp = grp.sort_values(["_tier", rank_col])

        names = (
            ";".join(grp[name_col].dropna().astype(str).unique())
            if name_col else ""
        )
        ids = (
            ";".join(grp[id_col].dropna().astype(str).unique())
            if id_col else ""
        )
        best_source = ""
        if source_col and not grp[source_col].isna().all():
            tier_vals = grp[source_col].dropna().map(GENE_SOURCE_RANK)
            if not tier_vals.empty:
                best_source = grp.loc[tier_vals.idxmin(), source_col]

        result[snp_id.lower()] = {
            "gene_names":  names,
            "gene_ids":    ids,
            "gene_source": best_source,
        }

    log.info("SNP→gene lookup built for %d rsIDs.", len(result))
    return result


# ---------------------------------------------------------------------------
# STAGE 2 – Merge gene annotations onto feature scores
# ---------------------------------------------------------------------------

def merge_genes_into_feature_scores(
    feature_scores_file: str,
    snp_gene_map: dict[str, dict],
) -> pd.DataFrame:
    """
    Annotate every row of *feature_scores_file* with gene information.

    All original rows and columns are preserved.  Non-SNP rows simply receive
    NaN in the added gene columns.

    Added columns
    -------------
    SNP1                : first (or only) rsID from 'feature' cell
    gene_names_SNP1     : semicolon-joined gene name(s) for SNP1 (best tier first)
    snp1_gene_ids       : semicolon-joined NCBI gene ID(s) for SNP1
    snp1_gene_source    : how genes were found — 'snp_db' | 'elink' | 'gene_db'
    SNP2                : second rsID (GxG rows only)
    gene_names_SNP2     : gene name(s) for SNP2
    snp2_gene_ids       : gene ID(s) for SNP2
    snp2_gene_source    : gene lookup tier for SNP2
    entrez_matched      : True if at least one rsID in this row had any gene hit
    """
    log.info("Loading feature scores file for annotation: %s", feature_scores_file)
    df = pd.read_csv(feature_scores_file)
    log.info("  %d rows × %d columns loaded.", *df.shape)

    if "feature" not in df.columns:
        raise ValueError("feature_scores_file must contain a 'feature' column.")

    def _annotate(feature: str) -> dict:
        rsids = _rsids_from_feature(feature)
        ftype = _feature_type(rsids)

        out = {
            "snp1"            : None,
            "snp1_gene_names" : None,
            "snp1_gene_ids"   : None,
            "snp1_gene_source": None,
            "snp2"            : None,
            "snp2_gene_names" : None,
            "snp2_gene_ids"   : None,
            "snp2_gene_source": None,
            "entrez_matched"  : False,
        }

        feat_lower = str(feature).strip().lower()
        if ftype == "other" and _is_hla_feature(feat_lower):
            hla_info = snp_gene_map.get(feat_lower, {})
            if hla_info:
                out["snp1"]             = _hla_feature_to_query(feat_lower)
                out["snp1_gene_names"]  = hla_info.get("gene_names")  or None
                out["snp1_gene_ids"]    = hla_info.get("gene_ids")    or None
                out["snp1_gene_source"] = hla_info.get("gene_source") or None
                out["entrez_matched"]   = True
            return out

        if ftype == "other":
            return out

        snp1      = rsids[0]
        snp1_info = snp_gene_map.get(snp1, {})
        out["snp1"]             = snp1
        out["snp1_gene_names"]  = snp1_info.get("gene_names")  or None
        out["snp1_gene_ids"]    = snp1_info.get("gene_ids")    or None
        out["snp1_gene_source"] = snp1_info.get("gene_source") or None
        if snp1_info:
            out["entrez_matched"] = True

        if ftype == "GxG" and len(rsids) >= 2:
            snp2      = rsids[1]
            snp2_info = snp_gene_map.get(snp2, {})
            out["snp2"]             = snp2
            out["snp2_gene_names"]  = snp2_info.get("gene_names")  or None
            out["snp2_gene_ids"]    = snp2_info.get("gene_ids")    or None
            out["snp2_gene_source"] = snp2_info.get("gene_source") or None
            if snp2_info:
                out["entrez_matched"] = True

        return out

    log.info("Annotating %d rows …", len(df))
    annotation_df = pd.DataFrame(
        df["feature"].map(_annotate).tolist(),
        index=df.index,
    )
    annotated = pd.concat([df, annotation_df], axis=1)

    annotated = annotated.rename(columns={
        "snp1":             "SNP1",
        "snp1_gene_names":  "gene_names_SNP1",
        "snp2":             "SNP2",
        "snp2_gene_names":  "gene_names_SNP2",
    })

    matched   = int(annotated["entrez_matched"].sum())
    snp_rows  = annotated["SNP1"].notna().sum()
    unmatched = int(snp_rows - matched)
    other     = int(annotated["SNP1"].isna().sum())

    source_counts = (
        annotated["snp1_gene_source"].value_counts().to_dict()
        if "snp1_gene_source" in annotated.columns else {}
    )
    log.info(
        "Annotation summary → rows with gene hit: %d "
        "(snp_db=%d elink=%d gene_db=%d hla=%d) | without hit: %d | non-queryable rows: %d",
        matched,
        source_counts.get("snp_db",      0),
        source_counts.get("elink",       0),
        source_counts.get("gene_db",     0),
        source_counts.get("hla_gene_db", 0),
        unmatched,
        other,
    )
    return annotated


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def run(
    filter_file: str,
    email: str,
    output_path: str,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], pd.DataFrame]:
    """
    Execute both pipeline stages and write all output files.

    Parameters
    ----------
    filter_file  : featureScoresReducedFinalModel.combined.filtered.csv.
                   Only features with non-zero 'coefs' are queried.
                   This file is also annotated in Stage 2.
    email        : NCBI Entrez email (required by NCBI policy).
    output_path  : Directory for all outputs (created if absent).
                   Existing gene_records.csv / unmatched_snps.csv are
                   loaded as a cache and extended with new results.

    Returns
    -------
    snp_records, gene_records, unmatched, feature_scores_annotated
    """
    Entrez.email = email
    out = Path(output_path)
    out.mkdir(parents=True, exist_ok=True)

    # ---- Stage 1A: load cache + extract queryable features -----------------
    cached_gene_records, cached_snp_records, already_queried = load_cache(out)

    snps, hla_map = extract_features_from_filter(filter_file)

    if not snps and not hla_map:
        log.warning("No queryable features (rsIDs or HLA alleles) found. Aborting.")
        return pd.DataFrame(), pd.DataFrame(), set(), pd.DataFrame()

    # ---- Stage 1B: query SNP db (three-tier) — skip cached SNPs -----------
    new_snp_records   = pd.DataFrame()
    new_gene_records_snp = pd.DataFrame()
    new_unmatched: set[str] = set()

    if snps:
        new_snp_records, new_gene_records_snp, new_unmatched = query_entrez_for_genes(
            snps, skip_snps=already_queried,
        )

    # ---- Stage 1B (HLA): query gene db for HLA allele features — skip cached
    new_gene_records_hla = query_hla_genes(hla_map, skip_keys=already_queried)

    # Merge new results with cache
    new_gene_records = pd.concat(
        [df for df in (new_gene_records_snp, new_gene_records_hla) if not df.empty],
        ignore_index=True,
    ) if (not new_gene_records_snp.empty or not new_gene_records_hla.empty) else pd.DataFrame()

    gene_records = pd.concat(
        [df for df in (cached_gene_records, new_gene_records) if not df.empty],
        ignore_index=True,
    ) if (not cached_gene_records.empty or not new_gene_records.empty) else pd.DataFrame()

    snp_records = pd.concat(
        [df for df in (cached_snp_records, new_snp_records) if not df.empty],
        ignore_index=True,
    ) if (not cached_snp_records.empty or not new_snp_records.empty) else pd.DataFrame()

    # Unmatched = previously unmatched (from cache) + newly unmatched,
    # minus any that were matched in a previous or current run
    cached_unmatched: set[str] = set()
    unmatched_p = out / "unmatched_snps.csv"
    if unmatched_p.exists():
        df_um = pd.read_csv(unmatched_p)
        if "rsid" in df_um.columns:
            cached_unmatched = set(df_um["rsid"].dropna().str.lower().unique())

    matched_rsids: set[str] = set()
    if not gene_records.empty and "snpId" in gene_records.columns:
        matched_rsids = set(gene_records["snpId"].dropna().str.lower().unique())

    unmatched = (cached_unmatched | new_unmatched) - matched_rsids

    # ---- Stage 1C: build lookup dict ---------------------------------------
    snp_gene_map = build_snp_gene_map(gene_records)

    # ---- Stage 1 outputs (write full merged files) -------------------------
    if not snp_records.empty:
        p = out / "g_features_snp_map.csv"
        snp_records.to_csv(p, index=False)
        log.info("SNP Entrez summary written to: %s", p)

    if not gene_records.empty:
        p = out / "gene_records.csv"
        gene_records.to_csv(p, index=False)
        log.info("Per-(feature × gene) records written to: %s", p)

    p = out / "unmatched_snps.csv"
    pd.DataFrame(sorted(unmatched), columns=["rsid"]).to_csv(p, index=False)
    log.info(
        "Unmatched rsIDs (%d) written to: %s",
        len(unmatched), p,
    )

    # ---- Stage 2: annotate filter_file -------------------------------------
    feature_scores_annotated = merge_genes_into_feature_scores(
        filter_file, snp_gene_map
    )

    p = out / "feature_scores_annotated.csv"
    feature_scores_annotated.to_csv(p, index=False)
    log.info("Annotated feature scores written to: %s", p)

    return snp_records, gene_records, unmatched, feature_scores_annotated


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Annotate a feature scores CSV with NCBI gene information.  "
            "Features are sourced from a penalized-model filter file (non-zero "
            "coefs only) and an optional supplementary SNP input file.  "
            "Previously queried SNPs are loaded from the output directory and "
            "skipped, enabling incremental re-runs."
        )
    )
    parser.add_argument(
        "--filter_file",
        type=str,
        required=True,
        help=(
            "featureScoresReducedFinalModel.combined.filtered.csv.  "
            "Must have 'coefs', 'feature', and (optionally) 'model' columns.  "
            "Only features with non-zero coefs are queried.  "
            "This file is annotated in Stage 2."
        ),
    )
    parser.add_argument(
        "--email",
        type=str,
        required=True,
        help="Your email address (required by NCBI Entrez usage policy).",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="./snpToGene_output",
        help=(
            "Output directory (created if absent). Default: ./snpToGene_output.  "
            "Existing gene_records.csv and unmatched_snps.csv are used as a "
            "cache; new results are appended."
        ),
    )

    args = parser.parse_args()

    log.info("filter_file  : %s", args.filter_file)
    log.info("output_path  : %s", args.output_path)

    run(
        filter_file  = args.filter_file,
        email        = args.email,
        output_path  = args.output_path,
    )
