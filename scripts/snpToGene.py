"""
snpToGene.py
------------
Two-stage pipeline driven entirely by the feature scores file — no
separate pre-filtered G/GxG file is required.

  STAGE 1 – Build feature→gene lookup
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Read the feature scores CSV and extract three classes of queryable features:

    rsID features   (G / GxG)   : 'rs12345' or 'rs12345,rs67890'
    HLA allele features          : 'DRB_5301' → queried as 'HLA-DRB*5301'
    Covariate features           : rows where model == 'covariate' — skipped
    Clinical / other features    : len(feature.split(',')) > 2, or no rsID
                                   and not HLA pattern — skipped

  Queries NCBI Entrez via three tiers (rsID) or gene-db direct search (HLA):

      g_features_snp_map.csv   – one row per rsID with full Entrez summary
      gene_records.csv         – one row per (feature × gene) pair
      unmatched_snps.csv       – rsIDs / HLA features with no gene hit

  STAGE 2 – Annotate feature scores
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Merge gene annotations back onto every row of the feature scores file.
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
    --feature_scores_file path/to/feature_scores.csv \\
    --email               you@institution.edu \\
    --output_path         ./snpToGene_output
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
# STAGE 1A – Extract queryable features from the feature scores file
# ---------------------------------------------------------------------------

def extract_features_from_scores(
    feature_scores_file: str,
) -> tuple[set[str], dict[str, str]]:
    """
    Extract queryable features from the full feature scores CSV.

    Classification rules (applied in order)
    ----------------------------------------
    1. Skip rows where ``model == 'covariate'``.  Covariate features
       (age, sex, etc.) appear in every model and should not be queried.
    2. Skip features where ``len(feature.split(',')) > 2``.  More than two
       comma-separated tokens signals a clinical feature encoding where
       index-0 is the trait name, not an rsID.
    3. Try to parse rsIDs via ``_rsids_from_feature``.  Valid hits → snps.
    4. If no rsIDs found and the feature matches the HLA allele pattern
       (exactly one ``_``, both parts alphanumeric, no commas) →
       convert to ``'HLA-<gene>*<allele>'`` and add to hla_map.
    5. All remaining features are clinical / non-genomic — skipped.

    Returns
    -------
    snps    : set[str]         lower-cased rsIDs to query via SNP db tiers
    hla_map : dict[str, str]   {original_feature_lower: 'HLA-DRB*5301'}
    """
    log.info("Loading feature scores file: %s", feature_scores_file)
    df = pd.read_csv(feature_scores_file)
    log.info("  %d rows × %d columns loaded.", *df.shape)

    if "feature" not in df.columns:
        raise ValueError(
            f"feature_scores_file must contain a 'feature' column. "
            f"Found columns: {df.columns.tolist()}"
        )

    # Identify covariate features — present in every model, never genomic
    covariate_features: set[str] = set()
    if "model" in df.columns:
        covariate_features = set(
            df.loc[df["model"] == "covariate", "feature"]
            .dropna().str.strip().str.lower()
        )
        if covariate_features:
            log.info("  Excluding %d covariate features.", len(covariate_features))

    snps:    set[str]        = set()
    hla_map: dict[str, str]  = {}
    n_clinical = 0

    for feat in df["feature"].dropna().unique():
        feat_str   = str(feat).strip()
        feat_lower = feat_str.lower()

        # Rule 1: skip covariates
        if feat_lower in covariate_features:
            continue

        # Rule 2: skip multi-token clinical features
        if len(feat_str.split(",")) > 2:
            n_clinical += 1
            continue

        # Rule 3: rsID(s)
        rsids = _rsids_from_feature(feat_str)
        if rsids:
            snps.update(rsids)
            continue

        # Rule 4: HLA allele notation
        if _is_hla_feature(feat_lower):
            hla_map[feat_lower] = _hla_feature_to_query(feat_lower)
            continue

        # Rule 5: everything else is clinical / non-genomic
        n_clinical += 1

    log.info(
        "  Extracted %d unique rsIDs, %d HLA features; %d clinical/other skipped.",
        len(snps), len(hla_map), n_clinical,
    )
    return snps, hla_map


# ---------------------------------------------------------------------------
# STAGE 1B – Query Entrez  (three-tier gene lookup)
# ---------------------------------------------------------------------------
#
# Tier 1 – SNP db GENES field  (most authoritative; gene_source = 'snp_db')
#   esearch(db="snp") → esummary → DocumentSummary.GENES
#   The SNP database directly annotates which genes a variant falls in or near.
#
# Tier 2 – elink SNP→gene  (gene_source = 'elink')
#   The SNP record exists but GENES is empty (common for intergenic / older
#   entries).  elink(dbfrom="snp", db="gene") traverses NCBI's cross-db
#   link table and usually recovers the neighbouring gene(s).
#
# Tier 3 – gene db direct search  (gene_source = 'gene_db')
#   The rsID is not in the SNP db at all (e.g. removed / merged rsIDs, or
#   variants only deposited in the gene db).  esearch(db="gene", term=rsid)
#   uses the gene db's full-text index which indexes associated rsIDs.
#   Results are returned in NCBI relevance order (most relevant gene first).
#
# gene_rank records the 1-based position of each gene within its tier's
# returned list, providing a lightweight proxy for relevance ordering.
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

    Parameters
    ----------
    queries    : ordered list of search terms to try
    feature_id : value to store in the snpId column of returned rows
    source     : gene_source label ('gene_db' or 'hla_gene_db')
    max_hits   : maximum number of gene UIDs to fetch

    Returns
    -------
    list of gene-row dicts (may be empty)
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

    Tries three progressively broader queries:
      1. ``'HLA-DRB*5301[gene name]'`` – field-tagged gene name search
      2. ``'HLA-DRB*5301'``            – free-text
      3. ``'HLA-DRB'``                 – gene symbol prefix (fallback)
    """
    gene_prefix = hla_query.split('*')[0]   # e.g. 'HLA-DRB'
    return _search_gene_db(
        queries=[f"{hla_query}[gene name]", hla_query, gene_prefix],
        feature_id=orig_feat,
        source="hla_gene_db",
    )


def query_hla_genes(hla_map: dict[str, str]) -> pd.DataFrame:
    """
    Query the NCBI gene database for each HLA allele feature.

    Parameters
    ----------
    hla_map : {original_feature_lower: 'HLA-DRB*5301', ...}
              Produced by extract_features_from_scores().

    Returns
    -------
    pd.DataFrame with columns: GENE_ID, NAME, snpId (=original feature key),
                                gene_source='hla_gene_db', gene_rank
    """
    if not hla_map:
        return pd.DataFrame()

    rows: list[dict] = []
    total = len(hla_map)
    for i, (orig_feat, hla_query) in enumerate(sorted(hla_map.items()), 1):
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
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """
    Query NCBI Entrez for gene associations for each rsID using a three-tier
    strategy that maximises recall while preserving source provenance.

    Tier 1 – SNP db GENES field  (gene_source = 'snp_db')
    Tier 2 – elink SNP→gene      (gene_source = 'elink')
    Tier 3 – gene db search      (gene_source = 'gene_db')

    Returns
    -------
    snp_records : pd.DataFrame
        Full SNP db DocumentSummary per rsID found in the SNP db
        ('snpId' and 'IDList' columns added).
    gene_records : pd.DataFrame
        One row per (rsID × gene).  Columns: GENE_ID, NAME, snpId,
        gene_source ('snp_db' | 'elink' | 'gene_db'), gene_rank.
    unmatched : set[str]
        rsIDs for which no gene was found via any tier.
    """
    snp_rows:  list[dict] = []
    gene_rows: list[dict] = []
    unmatched: set[str]   = set()

    total = len(snps)
    for i, snp in enumerate(sorted(snps), 1):
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
        # Tier 3: gene database direct search  (SNP not in SNP db at all,
        #         or elink also returned nothing)
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
        # Sort by tier rank then gene_rank so most relevant gene is first in
        # the joined string; fall back to original order when cols are absent.
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
        # Record the best (lowest-rank) source tier present for this SNP
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

    Note: ``feature_type`` (G / GxG / other) is used internally to drive
    annotation logic but is not written to the output because it is
    redundant with the existing ``model`` column.

    Parameters
    ----------
    feature_scores_file : str
        Full feature scores CSV (must have a 'feature' column).
    snp_gene_map : dict
        Lookup produced by build_snp_gene_map().

    Returns
    -------
    pd.DataFrame  –  annotated feature scores
    """
    log.info("Loading feature scores file: %s", feature_scores_file)
    df = pd.read_csv(feature_scores_file)
    log.info("  %d rows × %d columns loaded.", *df.shape)

    if "feature" not in df.columns:
        raise ValueError("feature_scores_file must contain a 'feature' column.")

    def _annotate(feature: str) -> dict:
        rsids = _rsids_from_feature(feature)
        ftype = _feature_type(rsids)   # used internally; not written to output

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

        # HLA allele features — looked up by original feature key (lowercased)
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

        # snp1 (always present for G and GxG)
        snp1      = rsids[0]
        snp1_info = snp_gene_map.get(snp1, {})
        out["snp1"]             = snp1
        out["snp1_gene_names"]  = snp1_info.get("gene_names")  or None
        out["snp1_gene_ids"]    = snp1_info.get("gene_ids")    or None
        out["snp1_gene_source"] = snp1_info.get("gene_source") or None
        if snp1_info:
            out["entrez_matched"] = True

        # snp2 (GxG only)
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

    # Rename to consistent output column names before any further references
    annotated = annotated.rename(columns={
        "snp1":             "SNP1",
        "snp1_gene_names":  "gene_names_SNP1",
        "snp2":             "SNP2",
        "snp2_gene_names":  "gene_names_SNP2",
    })

    matched   = int(annotated["entrez_matched"].sum())
    # Infer SNP vs other rows from the presence of a parsed SNP1 value;
    # feature_type is no longer written to the output.
    snp_rows  = annotated["SNP1"].notna().sum()
    unmatched = int(snp_rows - matched)
    other     = int(annotated["SNP1"].isna().sum())

    # Source-tier breakdown (genomic rows only)
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
    feature_scores_file: str,
    email: str,
    output_path: str,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], pd.DataFrame]:
    """
    Execute both pipeline stages and write all output files.

    Parameters
    ----------
    feature_scores_file : Full feature scores CSV.  rsIDs, HLA allele features,
                          and covariate exclusions are all derived from this file.
    email               : NCBI Entrez email (required by NCBI policy).
    output_path         : Directory for all outputs (created if absent).

    Returns
    -------
    snp_records, gene_records, unmatched, feature_scores_annotated
    """
    Entrez.email = email
    out = Path(output_path)
    out.mkdir(parents=True, exist_ok=True)

    # ---- Stage 1A: extract queryable features ------------------------------
    snps, hla_map = extract_features_from_scores(feature_scores_file)

    if not snps and not hla_map:
        log.warning("No queryable features (rsIDs or HLA alleles) found. Aborting.")
        return pd.DataFrame(), pd.DataFrame(), set(), pd.DataFrame()

    # ---- Stage 1B: query SNP db (three-tier) -------------------------------
    snp_records   = pd.DataFrame()
    gene_records_snp = pd.DataFrame()
    unmatched: set[str] = set()

    if snps:
        snp_records, gene_records_snp, unmatched = query_entrez_for_genes(snps)

    # ---- Stage 1B (HLA): query gene db for HLA allele features -------------
    gene_records_hla = query_hla_genes(hla_map)

    # Combine all gene records
    gene_records = pd.concat(
        [df for df in (gene_records_snp, gene_records_hla) if not df.empty],
        ignore_index=True,
    ) if (not gene_records_snp.empty or not gene_records_hla.empty) else pd.DataFrame()

    # ---- Stage 1C: build lookup dict ---------------------------------------
    snp_gene_map = build_snp_gene_map(gene_records)

    # ---- Stage 1 outputs ---------------------------------------------------
    if not snp_records.empty:
        p = out / "g_features_snp_map.csv"
        snp_records.to_csv(p, index=False)
        log.info("SNP Entrez summary written to: %s", p)

    if not gene_records.empty:
        p = out / "gene_records.csv"
        gene_records.to_csv(p, index=False)
        log.info("Per-(feature × gene) records written to: %s", p)

    # Write unmatched even when empty — downstream scripts can rely on its existence
    p = out / "unmatched_snps.csv"
    pd.DataFrame(sorted(unmatched), columns=["rsid"]).to_csv(p, index=False)
    log.info(
        "Unmatched rsIDs (%d) written to: %s  "
        "(ready for Ensembl / OpenTargets / dbSNP REST lookup)",
        len(unmatched), p,
    )

    # ---- Stage 2 -----------------------------------------------------------
    feature_scores_annotated = merge_genes_into_feature_scores(
        feature_scores_file, snp_gene_map
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
            "rsIDs, HLA allele features, and covariate exclusions are all "
            "derived automatically from the feature scores file."
        )
    )
    parser.add_argument(
        "--feature_scores_file",
        type=str,
        required=True,
        help=(
            "Full feature scores CSV.  Must have 'feature' and (optionally) "
            "'model' columns.  Covariate rows (model == 'covariate') are "
            "excluded from querying; all rows are preserved in the output."
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
        help="Output directory (created if absent). Default: ./snpToGene_output",
    )

    args = parser.parse_args()

    log.info("feature_scores_file : %s", args.feature_scores_file)
    log.info("output_path         : %s", args.output_path)

    run(
        feature_scores_file = args.feature_scores_file,
        email               = args.email,
        output_path         = args.output_path,
    )
