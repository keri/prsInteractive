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
# When POPULATION=separate the script loops automatically over both
# {cohort}/high_risk/ and {cohort}/low_control/ subdirectories, passing
# the matching --population_context to the literature enrichment script so:
#   high_risk   → standard risk-oriented queries
#   low_control → adds protective/prevention queries + unique protective columns
#
# Prerequisite: run cohort_bnmf_subtype_analysis.py first to produce
# feature_loadings.csv in the cluster output directory.
#
# ── Usage examples ────────────────────────────────────────────────────────────
#
# Single directory (auto-detect context from path):
#   CLUSTER_DIR=/path/to/cohortBnmf/cardio \
#   PHENOTYPES="type 2 diabetes" \
#   ./workflows/run_cluster_enrichment.sh
#
# Explicit high-risk context:
#   CLUSTER_DIR=/path/to/cohortBnmf/cardio/high_risk \
#   PHENOTYPES="type 2 diabetes" \
#   POPULATION_CONTEXT=high_risk \
#   ./workflows/run_cluster_enrichment.sh
#
# Separate populations (loops over both subdirs automatically):
#   BASE_COHORT_DIR=/path/to/cohortBnmf/cardio \
#   PHENOTYPES="type 2 diabetes" \
#   POPULATION=separate \
#   ./workflows/run_cluster_enrichment.sh
#
# ── Required environment variables ───────────────────────────────────────────
#   CLUSTER_DIR or BASE_COHORT_DIR — see below
#   PHENOTYPES    — comma-separated phenotype terms for literature queries
#
# ── Optional overrides (shown with defaults) ─────────────────────────────────
#   POPULATION          — 'single' (use CLUSTER_DIR) | 'separate' (loop over
#                         BASE_COHORT_DIR/{high_risk,low_control}); default: single
#   BASE_COHORT_DIR     — parent dir used when POPULATION=separate
#   POPULATION_CONTEXT  — 'auto'|'high_risk'|'low_control'|'both'; default: auto
#   TOP_N               — top features per cluster to query (default: 20)
#   TISSUES             — comma-separated GTEx tissue IDs (default: panel)
#   NCBI_API_KEY        — NCBI Entrez API key for 10 req/sec (default: 3/sec)
#   DISGENET_KEY        — DisGeNET API key (default: unauthenticated)
#   BACKGROUND_TABLE    — path to phenotype_background.csv from
#                         phenotype_background_analysis.py; when set, ALL terms
#                         in the table are used as phenotype query terms so every
#                         gene/SNP gets a score for each background term
#   NO_APPEND           — set to 1 to overwrite existing CSVs (default: append)
#   SKIP_GTEX           — set to 1 to skip GTEx enrichment (default: 0)
#   SKIP_LITERATURE     — set to 1 to skip literature enrichment (default: 0)

set -euo pipefail

# ── Resolve run mode ──────────────────────────────────────────────────────────
POPULATION="${POPULATION:-single}"

if [[ "${POPULATION}" == "separate" ]]; then
  : "${BASE_COHORT_DIR:?ERROR: BASE_COHORT_DIR must be set when POPULATION=separate}"
  DIRS=("${BASE_COHORT_DIR}/high_risk" "${BASE_COHORT_DIR}/low_control")
  CONTEXTS=("high_risk" "low_control")
else
  : "${CLUSTER_DIR:?ERROR: CLUSTER_DIR env var must be set}"
  DIRS=("${CLUSTER_DIR}")
  CONTEXTS=("${POPULATION_CONTEXT:-auto}")
fi

: "${PHENOTYPES:?ERROR: PHENOTYPES env var must be set}"

# ── Optional with defaults ────────────────────────────────────────────────────
TOP_N="${TOP_N:-20}"
NCBI_API_KEY="${NCBI_API_KEY:-}"
DISGENET_KEY="${DISGENET_KEY:-}"
BACKGROUND_TABLE="${BACKGROUND_TABLE:-}"
SKIP_GTEX="${SKIP_GTEX:-0}"
SKIP_LITERATURE="${SKIP_LITERATURE:-0}"
NO_APPEND="${NO_APPEND:-0}"

SCRIPT_DIR="$(cd "$(dirname "$0")/../scripts" && pwd)"

echo "============================================================"
echo "CLUSTER ENRICHMENT PIPELINE"
echo "============================================================"
echo "  POPULATION    : ${POPULATION}"
echo "  PHENOTYPES    : ${PHENOTYPES}"
echo "  TOP_N         : ${TOP_N}"
echo "  Directories   : ${DIRS[*]}"
echo ""

# ── Build optional argument lists ─────────────────────────────────────────────
GTEX_ARGS=()
LIT_ARGS_BASE=()

[[ -n "${NCBI_API_KEY}" ]] && GTEX_ARGS+=(--ncbi_api_key "${NCBI_API_KEY}")
[[ -n "${NCBI_API_KEY}" ]]    && LIT_ARGS_BASE+=(--ncbi_api_key  "${NCBI_API_KEY}")
[[ -n "${DISGENET_KEY}" ]]    && LIT_ARGS_BASE+=(--disgenet_api_key "${DISGENET_KEY}")
[[ -n "${TISSUES:-}" ]]       && GTEX_ARGS+=(--tissues "${TISSUES}")
[[ -n "${BACKGROUND_TABLE}" ]] && LIT_ARGS_BASE+=(--background_table "${BACKGROUND_TABLE}")
[[ "${NO_APPEND}" == "1" ]]   && LIT_ARGS_BASE+=(--no_append)

# ── Per-directory loop ────────────────────────────────────────────────────────
for i in "${!DIRS[@]}"; do
  DIR="${DIRS[$i]}"
  CTX="${CONTEXTS[$i]}"

  echo "============================================================"
  echo "  Directory : ${DIR}"
  echo "  Context   : ${CTX}"
  echo "============================================================"

  if [[ ! -f "${DIR}/feature_loadings.csv" ]]; then
    echo "  [SKIP] feature_loadings.csv not found in ${DIR} — skipping"
    echo ""
    continue
  fi

  # ── Step 1: GTEx enrichment ────────────────────────────────────────────
  if [[ "${SKIP_GTEX}" != "1" ]]; then
    echo "------------------------------------------------------------"
    echo "Step 1: GTEx tissue expression + eQTL enrichment"
    echo "------------------------------------------------------------"
    python "${SCRIPT_DIR}/enrich_clusters_gtex.py" \
      --cluster_dir "${DIR}" \
      --top_n       "${TOP_N}" \
      "${GTEX_ARGS[@]+"${GTEX_ARGS[@]}"}"
    echo ""
  else
    echo "  [SKIP] GTEx enrichment (SKIP_GTEX=1)"
  fi

  # ── Step 2: Literature enrichment ─────────────────────────────────────
  if [[ "${SKIP_LITERATURE}" != "1" ]]; then
    echo "------------------------------------------------------------"
    echo "Step 2: Literature enrichment"
    echo "  DisGeNET + Open Targets + PubTator3 + PubMed + ClinVar"
    if [[ "${CTX}" == "low_control" ]]; then
      echo "  + Protective/prevention queries (low_control context)"
    fi
    echo "------------------------------------------------------------"
    python "${SCRIPT_DIR}/enrich_clusters_literature.py" \
      --cluster_dir        "${DIR}" \
      --phenotypes         "${PHENOTYPES}" \
      --top_n              "${TOP_N}" \
      --population_context "${CTX}" \
      "${LIT_ARGS_BASE[@]+"${LIT_ARGS_BASE[@]}"}"
    echo ""
  else
    echo "  [SKIP] Literature enrichment (SKIP_LITERATURE=1)"
  fi

  # ── Per-dir summary ────────────────────────────────────────────────────
  echo "  Outputs:"
  [[ "${SKIP_GTEX}"       != "1" ]] && echo "    GTEx       : ${DIR}/enrichment/gtex/"
  [[ "${SKIP_LITERATURE}" != "1" ]] && echo "    Literature : ${DIR}/enrichment/literature/"
  echo ""
done

# ── Final summary ─────────────────────────────────────────────────────────────
echo "============================================================"
echo "Enrichment Output File Reference"
echo "============================================================"
if [[ "${SKIP_GTEX}" != "1" ]]; then
  echo "  GTEx (per directory):"
  echo "    gtex_gene_expression.csv         — gene × tissue median TPM"
  echo "    gtex_eqtl_evidence.csv           — SNP → gene eQTL effect sizes per tissue"
  echo "    gtex_cluster_tissue_summary.csv  — cluster × tissue enrichment scores"
  echo "    gtex_cluster_tissue_heatmap.png  — heatmap of top tissues per cluster"
  echo ""
fi
if [[ "${SKIP_LITERATURE}" != "1" ]]; then
  echo "  Literature (per directory):"
  echo "    gene_evidence_table.csv          — term × phenotype evidence scores"
  echo "      risk columns  : disgenet_score, open_targets_score, pubtator3_count,"
  echo "                       pubmed_count, clinvar_sig, combined_score"
  echo "      prot columns  : protective_count, protective_sentiment,"
  echo "                       unique_protective, protective_titles"
  echo "                       (present only when context = low_control)"
  echo "    cluster_evidence_ranked.csv      — features ranked by evidence per cluster"
  echo "    cluster_literature_summary.csv   — dominant phenotype per cluster"
  echo ""
fi
echo "Done."
