#!/usr/bin/env bash
# workflows/run_cluster_enrichment.sh
#
# Post-NMF cluster enrichment pipeline.
#
# Runs two enrichment scripts against a bNMF output directory:
#   1. enrich_clusters_gtex.py      — GTEx tissue expression + eQTL evidence
#   2. enrich_clusters_literature.py — DisGeNET + Open Targets + PubTator3 +
#                                      PubMed + ClinVar
#
# Open Targets and PubTator3 require no API key and are queried automatically.
# DisGeNET and NCBI run unauthenticated by default (lower rate limits).
#
# Prerequisite: run cohort_bnmf_subtype_analysis.py first to produce
# feature_loadings.csv in the cluster output directory.
#
# Usage (per-cohort):
#   CLUSTER_DIR=/path/to/cohortBnmf/cardio \
#   PHENOTYPES="type 2 diabetes,epilepsy,cardio" \
#   ./workflows/run_cluster_enrichment.sh
#
# Usage (combined):
#   CLUSTER_DIR=/path/to/combinedCohortBnmf \
#   PHENOTYPES="type 2 diabetes,epilepsy,cardio" \
#   ./workflows/run_cluster_enrichment.sh
#
# ── Required environment variables ───────────────────────────────────────────
#   CLUSTER_DIR   — path to bNMF output directory containing feature_loadings.csv
#   PHENOTYPES    — comma-separated phenotype terms for literature queries
#
# ── Optional overrides (shown with defaults) ─────────────────────────────────
#   TOP_N             — top features per cluster to query (default: 20)
#   TISSUES           — comma-separated GTEx tissue IDs (default: cardiometabolic panel)
#   NCBI_API_KEY      — NCBI Entrez API key for 10 req/sec (default: none = 3 req/sec)
#   DISGENET_KEY      — DisGeNET API key (default: none = unauthenticated)
#   SKIP_GTEX         — set to 1 to skip GTEx enrichment (default: 0)
#   SKIP_LITERATURE   — set to 1 to skip literature enrichment (default: 0)

set -euo pipefail

# ── Required ──────────────────────────────────────────────────────────────────
: "${CLUSTER_DIR:?ERROR: CLUSTER_DIR env var must be set}"
: "${PHENOTYPES:?ERROR: PHENOTYPES env var must be set}"

# ── Optional with defaults ────────────────────────────────────────────────────
TOP_N="${TOP_N:-20}"
NCBI_API_KEY="${NCBI_API_KEY:-}"
DISGENET_KEY="${DISGENET_KEY:-}"
SKIP_GTEX="${SKIP_GTEX:-0}"
SKIP_LITERATURE="${SKIP_LITERATURE:-0}"

SCRIPT_DIR="$(cd "$(dirname "$0")/../scripts" && pwd)"

echo "============================================================"
echo "CLUSTER ENRICHMENT PIPELINE"
echo "============================================================"
echo "  CLUSTER_DIR  : ${CLUSTER_DIR}"
echo "  PHENOTYPES   : ${PHENOTYPES}"
echo "  TOP_N        : ${TOP_N}"
echo ""

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${CLUSTER_DIR}/feature_loadings.csv" ]]; then
  echo "  ERROR: feature_loadings.csv not found in ${CLUSTER_DIR}"
  echo "         Run cohort_bnmf_subtype_analysis.py --mode combined first."
  exit 1
fi

# ── Build optional argument lists ─────────────────────────────────────────────
GTEX_ARGS=()
LIT_ARGS=()

if [[ -n "${NCBI_API_KEY}" ]]; then
  GTEX_ARGS+=(--ncbi_api_key "${NCBI_API_KEY}")
  LIT_ARGS+=(--ncbi_api_key  "${NCBI_API_KEY}")
fi

if [[ -n "${DISGENET_KEY}" ]]; then
  LIT_ARGS+=(--disgenet_api_key "${DISGENET_KEY}")
fi

if [[ -n "${TISSUES:-}" ]]; then
  GTEX_ARGS+=(--tissues "${TISSUES}")
fi

# ── Step 1: GTEx enrichment ───────────────────────────────────────────────────
if [[ "${SKIP_GTEX}" != "1" ]]; then
  echo "------------------------------------------------------------"
  echo "Step 1: GTEx tissue expression + eQTL enrichment"
  echo "------------------------------------------------------------"
  python "${SCRIPT_DIR}/enrich_clusters_gtex.py" \
    --cluster_dir "${CLUSTER_DIR}" \
    --top_n       "${TOP_N}" \
    "${GTEX_ARGS[@]+"${GTEX_ARGS[@]}"}"
  echo ""
else
  echo "  [SKIP] GTEx enrichment (SKIP_GTEX=1)"
fi

# ── Step 2: Literature enrichment ─────────────────────────────────────────────
if [[ "${SKIP_LITERATURE}" != "1" ]]; then
  echo "------------------------------------------------------------"
  echo "Step 2: Literature enrichment (DisGeNET + Open Targets + PubTator3 + PubMed + ClinVar)"
  echo "------------------------------------------------------------"
  python "${SCRIPT_DIR}/enrich_clusters_literature.py" \
    --cluster_dir "${CLUSTER_DIR}" \
    --phenotypes  "${PHENOTYPES}" \
    --top_n       "${TOP_N}" \
    "${LIT_ARGS[@]+"${LIT_ARGS[@]}"}"
  echo ""
else
  echo "  [SKIP] Literature enrichment (SKIP_LITERATURE=1)"
fi

# ── Summary ────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Enrichment Outputs"
echo "============================================================"
if [[ "${SKIP_GTEX}" != "1" ]]; then
  echo "  GTEx results      : ${CLUSTER_DIR}/enrichment/gtex/"
  echo "    gtex_gene_expression.csv      — gene × tissue median TPM"
  echo "    gtex_eqtl_evidence.csv        — SNP → gene eQTL effect sizes per tissue"
  echo "    gtex_cluster_tissue_summary.csv — cluster × tissue enrichment scores"
  echo "    gtex_cluster_tissue_heatmap.png — heatmap of top tissues per cluster"
  echo ""
fi
if [[ "${SKIP_LITERATURE}" != "1" ]]; then
  echo "  Literature results: ${CLUSTER_DIR}/enrichment/literature/"
  echo "    gene_evidence_table.csv        — all term × phenotype evidence scores"
  echo "      columns: disgenet_score, open_targets_score, pubtator3_count,"
  echo "               pubmed_count, clinvar_sig, combined_score"
  echo "    cluster_evidence_ranked.csv    — per-cluster features ranked by evidence"
  echo "    cluster_literature_summary.csv — dominant phenotype per cluster"
  echo ""
fi
echo "Done."
