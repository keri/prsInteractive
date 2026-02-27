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
PARAM_EUTILS = {"usehistory": "Y", "retmode": "xml", "retmax": 10}

_RS_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


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
# STAGE 1B – Query Entrez
# ---------------------------------------------------------------------------

def _get_id_list(snp_id: str) -> list[str]:
    """Return the NCBI UID list for *snp_id*."""
    handle = Entrez.esearch(db=DB, term=snp_id, **PARAM_EUTILS)
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])


def _parse_esummary(dict_obj: dict) -> dict:
    """Extract the first DocumentSummary; return {} on failure."""
    try:
        return dict_obj["DocumentSummarySet"]["DocumentSummary"][0]
    except (KeyError, IndexError) as exc:
        log.warning("Could not parse eSummary: %s", exc)
        return {}


def query_entrez_for_genes(
    snps: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """
    Query NCBI Entrez SNP database for each rsID.

    Returns
    -------
    snp_records : pd.DataFrame
        Full Entrez DocumentSummary per matched rsID ('snpId' and 'IDList' added).
    gene_records : pd.DataFrame
        One row per (rsID × gene).  Columns include GENE_ID, NAME, ACC, snpId.
    unmatched : set[str]
        rsIDs with no Entrez UID, failed requests, or empty GENES list.
    """
    snp_rows:  list[dict] = []
    gene_rows: list[dict] = []
    unmatched: set[str]   = set()

    total = len(snps)
    for i, snp in enumerate(sorted(snps), 1):
        log.info("[%d/%d] Querying Entrez for %s …", i, total, snp)

        # --- esearch --------------------------------------------------------
        try:
            id_list = _get_id_list(snp)
        except Exception as exc:
            log.error("  esearch failed: %s", exc)
            unmatched.add(snp)
            continue

        if not id_list:
            log.debug("  No UIDs returned.")
            unmatched.add(snp)
            continue

        # --- esummary -------------------------------------------------------
        try:
            handle       = Entrez.esummary(db=DB, id=id_list[0], **PARAM_EUTILS)
            summary_dict = Entrez.read(handle)
            handle.close()
            record = _parse_esummary(summary_dict)
        except Exception as exc:
            log.error("  esummary failed (uid=%s): %s", id_list[0], exc)
            unmatched.add(snp)
            continue

        if not record:
            unmatched.add(snp)
            continue

        # --- SNP row --------------------------------------------------------
        row           = dict(record)
        row["snpId"]  = snp
        row["IDList"] = id_list
        snp_rows.append(row)

        # --- Gene rows (one per gene) ----------------------------------------
        genes = record.get("GENES", [])
        if genes:
            for g in genes:
                gene_row          = dict(g)
                gene_row["snpId"] = snp
                gene_rows.append(gene_row)
        else:
            log.debug("  Record found but GENES list is empty — added to unmatched.")
            unmatched.add(snp)

    snp_records  = pd.DataFrame(snp_rows)
    gene_records = pd.DataFrame(gene_rows)

    log.info(
        "Entrez complete → matched: %d SNP records | unmatched: %d",
        len(snp_rows), len(unmatched),
    )
    return snp_records, gene_records, unmatched


# ---------------------------------------------------------------------------
# STAGE 1C – Build the rsID → gene lookup dict
# ---------------------------------------------------------------------------

def build_snp_gene_map(gene_records: pd.DataFrame) -> dict[str, dict]:
    """
    Collapse the per-(rsID × gene) DataFrame into a dict keyed by lower-case rsID.

    When a SNP maps to multiple genes, names and IDs are joined with ';'.

    Returns
    -------
    {
      'rs1234': {'gene_names': 'BRCA1',       'gene_ids': '672'},
      'rs5678': {'gene_names': 'TP53;MDM2',   'gene_ids': '7157;4193'},
      ...
    }
    """
    if gene_records.empty:
        log.warning("gene_records is empty — snp_gene_map will be empty.")
        return {}

    # Normalise column names (Entrez returns mixed case; upper for safety)
    gr = gene_records.copy()
    gr.columns = [c.upper() for c in gr.columns]

    # After upper-casing, 'snpId' → 'SNPID'
    name_col = "NAME"    if "NAME"    in gr.columns else None
    id_col   = "GENE_ID" if "GENE_ID" in gr.columns else None

    result: dict[str, dict] = {}
    for snp_id, grp in gr.groupby("SNPID"):
        names = (
            ";".join(grp[name_col].dropna().astype(str).unique())
            if name_col else ""
        )
        ids = (
            ";".join(grp[id_col].dropna().astype(str).unique())
            if id_col else ""
        )
        result[snp_id.lower()] = {"gene_names": names, "gene_ids": ids}

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
    feature_type    : 'G' | 'GxG' | 'other'
    snp1            : first (or only) rsID from 'feature' cell
    snp1_gene_names : semicolon-joined gene name(s) for snp1
    snp1_gene_ids   : semicolon-joined NCBI gene ID(s) for snp1
    snp2            : second rsID (GxG rows only)
    snp2_gene_names : gene name(s) for snp2
    snp2_gene_ids   : gene ID(s) for snp2
    entrez_matched  : True if at least one rsID in this row had an Entrez hit

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
        ftype = _feature_type(rsids)

        out = {
            "feature_type"   : ftype,
            "snp1"           : None,
            "snp1_gene_names": None,
            "snp1_gene_ids"  : None,
            "snp2"           : None,
            "snp2_gene_names": None,
            "snp2_gene_ids"  : None,
            "entrez_matched" : False,
        }

        if ftype == "other":
            return out

        # snp1 (always present for G and GxG)
        snp1      = rsids[0]
        snp1_info = snp_gene_map.get(snp1, {})
        out["snp1"]            = snp1
        out["snp1_gene_names"] = snp1_info.get("gene_names") or None
        out["snp1_gene_ids"]   = snp1_info.get("gene_ids")   or None
        if snp1_info:
            out["entrez_matched"] = True

        # snp2 (GxG only)
        if ftype == "GxG" and len(rsids) >= 2:
            snp2      = rsids[1]
            snp2_info = snp_gene_map.get(snp2, {})
            out["snp2"]            = snp2
            out["snp2_gene_names"] = snp2_info.get("gene_names") or None
            out["snp2_gene_ids"]   = snp2_info.get("gene_ids")   or None
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
    snp_rows  = annotated["feature_type"].isin(["G", "GxG"]).sum()
    unmatched = int(snp_rows - matched)
    other     = int((annotated["feature_type"] == "other").sum())

    log.info(
        "Annotation summary → G/GxG rows with gene hit: %d | without hit: %d | other rows: %d",
        matched, unmatched, other,
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
