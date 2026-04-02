"""
snpToGene.py
------------
Two-stage pipeline:

  STAGE 1 – Build SNP→gene lookup
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Read a pre-filtered G / GxG features file (no HLA, no covariates, no
  E-feature labels — just bare rsIDs and comma-separated rsID pairs).
  Query NCBI Entrez for each unique rsID and build a lookup table:

      g_features_snp_map.csv   – one row per rsID with full Entrez summary
      gene_records.csv         – one row per (rsID × gene) pair
      unmatched_snps.csv       – rsIDs with no Entrez hit (for downstream
                                 searches: Ensembl, OpenTargets, dbSNP REST…)

  STAGE 2 – Annotate feature scores
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Read the full feature_scores file (may contain any model / feature type)
  and merge gene annotations from Stage 1 back onto every row via the
  'feature' column.  Non-SNP rows (covariates, HLA, E-labels) receive NaN
  in all gene columns so the full dataset is preserved.

  Output:
      feature_scores_annotated.csv – original feature_scores with added columns:
        feature_type        : 'G' | 'GxG' | 'other'
        snp1                : first (or only) rsID parsed from feature
        snp1_gene_names     : semicolon-joined gene names for snp1
        snp1_gene_ids       : semicolon-joined NCBI gene IDs for snp1
        snp2                : second rsID (GxG only, else NaN)
        snp2_gene_names     : semicolon-joined gene names for snp2
        snp2_gene_ids       : semicolon-joined NCBI gene IDs for snp2
        entrez_matched      : True / False (any rsID in feature had a hit)

Usage
-----
python snpToGene.py \\
    --g_features_file     path/to/g_and_gxg_features.csv \\
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
GENE_SOURCE_RANK = {"snp_db": 1, "elink": 2, "gene_db": 3}


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


# ---------------------------------------------------------------------------
# STAGE 1A – Parse the clean G / GxG features file
# ---------------------------------------------------------------------------

def extract_snps_from_g_features(g_features_file: str) -> set[str]:
    """
    Load the pre-filtered G / GxG features CSV and return every unique rsID.

    The file is trusted to contain only SNP features — no filtering is applied.
    Both single rsIDs and comma-separated epistatic pairs are handled.

    Parameters
    ----------
    g_features_file : str
        CSV with at minimum a 'feature' column.

    Returns
    -------
    set[str]  –  lower-cased rsIDs, e.g. {'rs1234', 'rs5678'}
    """
    log.info("Loading G/GxG features file: %s", g_features_file)
    df = pd.read_csv(g_features_file)
    log.info("  %d rows × %d columns loaded.", *df.shape)

    if "feature" not in df.columns:
        raise ValueError(
            f"G features file must contain a 'feature' column. "
            f"Found columns: {df.columns.tolist()}"
        )

    snps: set[str] = set()
    for feat in df["feature"].dropna().unique():
        snps.update(_rsids_from_feature(feat))

    log.info("  %d unique rsIDs extracted.", len(snps))
    return snps


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


def _search_gene_db_for_snp(snp_id: str) -> list[dict]:
    """
    Search the gene database directly for *snp_id* (Tier 3 fallback).

    Tries two progressively broader queries:
      1. ``<rsid>[SNP]``    – field-tagged rsID search
      2. ``<rsid>``         – free-text (catches gene entries that mention the rsID)

    Returns up to 5 gene rows ranked by NCBI relevance (gene_rank 1 = best).
    """
    for term in (f"{snp_id}[SNP]", snp_id):
        try:
            handle  = Entrez.esearch(db=DB_GENE, term=term, **PARAM_GENE)
            result  = Entrez.read(handle)
            handle.close()
            id_list = result.get("IdList", [])
        except Exception as exc:
            log.warning("  gene db esearch (%r) failed: %s", term, exc)
            id_list = []
        if id_list:
            log.debug("  gene db search %r → %d hits", term, len(id_list))
            return _fetch_gene_summaries(id_list[:5], snp_id, source="gene_db")
    return []


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
    snp1                : first (or only) rsID from 'feature' cell
    snp1_gene_names     : semicolon-joined gene name(s) for snp1 (best tier first)
    snp1_gene_ids       : semicolon-joined NCBI gene ID(s) for snp1
    snp1_gene_source    : how genes were found — 'snp_db' | 'elink' | 'gene_db'
    snp2                : second rsID (GxG rows only)
    snp2_gene_names     : gene name(s) for snp2
    snp2_gene_ids       : gene ID(s) for snp2
    snp2_gene_source    : gene lookup tier for snp2
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

    matched   = int(annotated["entrez_matched"].sum())
    # Infer SNP vs other rows from the presence of a parsed snp1 value;
    # feature_type is no longer written to the output.
    snp_rows  = annotated["snp1"].notna().sum()
    unmatched = int(snp_rows - matched)
    other     = int(annotated["snp1"].isna().sum())

    # Source-tier breakdown (genomic rows only)
    source_counts = (
        annotated["snp1_gene_source"].value_counts().to_dict()
        if "snp1_gene_source" in annotated.columns else {}
    )
    log.info(
        "Annotation summary → SNP rows with gene hit: %d "
        "(snp_db=%d elink=%d gene_db=%d) | without hit: %d | non-SNP rows: %d",
        matched,
        source_counts.get("snp_db",  0),
        source_counts.get("elink",   0),
        source_counts.get("gene_db", 0),
        unmatched,
        other,
    )
    return annotated


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def run(
    g_features_file: str,
    feature_scores_file: str,
    email: str,
    output_path: str,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], pd.DataFrame]:
    """
    Execute both pipeline stages and write all output files.

    Parameters
    ----------
    g_features_file     : CSV with only G / GxG features (pre-filtered).
    feature_scores_file : Full feature scores CSV to annotate.
    email               : NCBI Entrez email.
    output_path         : Directory for all outputs (created if absent).

    Returns
    -------
    snp_records, gene_records, unmatched, feature_scores_annotated
    """
    Entrez.email = email
    out = Path(output_path)
    out.mkdir(parents=True, exist_ok=True)

    # ---- Stage 1A ----------------------------------------------------------
    snps = extract_snps_from_g_features(g_features_file)

    if not snps:
        log.warning("No valid rsIDs found in G features file. Aborting.")
        return pd.DataFrame(), pd.DataFrame(), set(), pd.DataFrame()

    # ---- Stage 1B ----------------------------------------------------------
    snp_records, gene_records, unmatched = query_entrez_for_genes(snps)

    # ---- Stage 1C ----------------------------------------------------------
    snp_gene_map = build_snp_gene_map(gene_records)

    # ---- Stage 1 outputs ---------------------------------------------------
    if not snp_records.empty:
        p = out / "g_features_snp_map.csv"
        snp_records.to_csv(p, index=False)
        log.info("SNP Entrez summary written to:   %s", p)

    if not gene_records.empty:
        p = out / "gene_records.csv"
        gene_records.to_csv(p, index=False)
        log.info("Per-(rsID × gene) records written to: %s", p)

    # Write unmatched even when empty — downstream scripts can rely on its existence
    p = out / "unmatched_snps.csv"
    pd.DataFrame(sorted(unmatched), columns=["rsid"]).to_csv(p, index=False)
    log.info(
        "Unmatched SNPs (%d) written to: %s  "
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
            "Build an SNP→gene lookup from a clean G/GxG features file, "
            "then annotate a full feature scores CSV with gene information."
        )
    )
    parser.add_argument(
        "--g_features_file",
        type=str,
        required=True,
        help=(
            "CSV containing only G and GxG features (pre-filtered: no HLA, "
            "no covariates, no E-feature labels).  Must have a 'feature' column."
        ),
    )
    parser.add_argument(
        "--feature_scores_file",
        type=str,
        required=True,
        help=(
            "Full feature scores CSV.  All rows are preserved; gene columns "
            "are NaN for non-SNP features.  Must have a 'feature' column."
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

    log.info("g_features_file     : %s", args.g_features_file)
    log.info("feature_scores_file : %s", args.feature_scores_file)
    log.info("output_path         : %s", args.output_path)

    run(
        g_features_file     = args.g_features_file,
        feature_scores_file = args.feature_scores_file,
        email               = args.email,
        output_path         = args.output_path,
    )
