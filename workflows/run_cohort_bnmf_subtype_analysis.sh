#!/usr/bin/env bash
# workflows/run_cohort_bnmf_subtype_analysis.sh
#
# Run bNMF subtype analysis using outputs from the cohort analysis pipeline
# (clinical + genomic features) applied to high-risk cases and low-control
# controls in each cohort.
#
# Two modes are available (controlled by MODE env var, default: both):
#   per_cohort — independent bNMF per cohort (clusters not cross-comparable)
#   combined   — single bNMF on union of all cohorts with cohort-prefixed features,
#                producing globally comparable clusters + pathway analysis outputs
#   both       — run per_cohort then combined (default)
#
# Prerequisite: run_cohort_analysis_pipeline.sh must have been run first to
# produce the genomic and clinical feature files consumed by this workflow.
#
# Usage:
#   PHENO_DATA=/path/to/results \
#   RAW_FEATURES_FILE=/path/to/participant_environment.csv \
#   ./workflows/run_cohort_bnmf_subtype_analysis.sh
#
# ── Required environment variables ───────────────────────────────────────────
#   PHENO_DATA           — path to phenotype results directory
#   RAW_FEATURES_FILE    — path to raw (untransformed) clinical features CSV
#
# ── Optional overrides (shown with defaults) ─────────────────────────────────
#   MODE                  — analysis mode: per_cohort|combined|both (default: both)
#   FILTER_STRATEGY       — genomic filter strategy: tiered|strict|majority|differential
#                           default: tiered
#   USE_SET               — individual set: holdout|validation|all
#                           default: holdout
#   K_MIN                 — minimum k for NMF consensus       (default: 2)
#   K_MAX                 — maximum k for NMF consensus       (default: 6)
#   N_RUNS                — NMF restarts per k                (default: 30)
#   K_SELECT_METHOD       — optimal-k method: elbow|max       (default: elbow)
#   L1_RATIO              — L1 sparsity ratio                 (default: 0.1)
#   ALPHA_W               — alpha_W regularisation            (default: 0.1)
#   MIN_EFFECT_SIZE       — min |rank-biserial r| for clinical (default: 0.0)
#   SPECIFICITY_TIER_MAX  — max genomic tier (1-3)            (default: 3)
#   HIGH_RISK_THRESHOLD   — bin > this = high-risk            (default: 8)
#   LOW_CONTROL_THRESHOLD — bin < this = low-control          (default: 3)
#   COHORTS               — space-separated cohort names to run (default: all)
#   NO_HOLDOUT_PRS        — set to 1 to exclude holdout PRS files (default: 0)

set -euo pipefail

# ── Required ──────────────────────────────────────────────────────────────────
: "${PHENO_DATA:?ERROR: PHENO_DATA env var must be set}"
: "${RAW_FEATURES_FILE:?ERROR: RAW_FEATURES_FILE env var must be set}"

# ── Optional with defaults ────────────────────────────────────────────────────
MODE="${MODE:-both}"
FILTER_STRATEGY="${FILTER_STRATEGY:-tiered}"
USE_SET="${USE_SET:-holdout}"
K_MIN="${K_MIN:-2}"
K_MAX="${K_MAX:-6}"
N_RUNS="${N_RUNS:-30}"
K_SELECT_METHOD="${K_SELECT_METHOD:-elbow}"
L1_RATIO="${L1_RATIO:-0.1}"
ALPHA_W="${ALPHA_W:-0.1}"
MIN_EFFECT_SIZE="${MIN_EFFECT_SIZE:-0.0}"
SPECIFICITY_TIER_MAX="${SPECIFICITY_TIER_MAX:-3}"
HIGH_RISK_THRESHOLD="${HIGH_RISK_THRESHOLD:-8}"
LOW_CONTROL_THRESHOLD="${LOW_CONTROL_THRESHOLD:-3}"
NO_HOLDOUT_PRS="${NO_HOLDOUT_PRS:-0}"

SCORES="${PHENO_DATA}/scores"

echo "============================================================"
echo "COHORT bNMF SUBTYPE ANALYSIS"
echo "============================================================"
echo "  PHENO_DATA            : ${PHENO_DATA}"
echo "  RAW_FEATURES_FILE     : ${RAW_FEATURES_FILE}"
echo "  MODE                  : ${MODE}"
echo "  FILTER_STRATEGY       : ${FILTER_STRATEGY}"
echo "  USE_SET               : ${USE_SET}"
echo "  K_MIN / K_MAX         : ${K_MIN} / ${K_MAX}"
echo "  N_RUNS                : ${N_RUNS}"
echo "  K_SELECT_METHOD       : ${K_SELECT_METHOD}"
echo "  L1_RATIO / ALPHA_W    : ${L1_RATIO} / ${ALPHA_W}"
echo "  HIGH_RISK_THRESHOLD   : bin > ${HIGH_RISK_THRESHOLD}"
echo "  LOW_CONTROL_THRESHOLD : bin < ${LOW_CONTROL_THRESHOLD}"
echo "  MIN_EFFECT_SIZE       : ${MIN_EFFECT_SIZE}"
echo "  SPECIFICITY_TIER_MAX  : ${SPECIFICITY_TIER_MAX}"
if [[ -n "${COHORTS:-}" ]]; then
  echo "  COHORTS               : ${COHORTS}"
else
  echo "  COHORTS               : all"
fi
echo ""

# ── Validate prerequisite files ───────────────────────────────────────────────
echo "Checking prerequisite files..."

GENOMIC_FILE="${SCORES}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.Filtered.${FILTER_STRATEGY}.csv"
CLINICAL_FILE="${SCORES}/cohortClinicalComparison/clinical_comparison_all_${USE_SET}.csv"

if [[ ! -f "${GENOMIC_FILE}" ]]; then
  echo "  [WARN] Genomic features file not found:"
  echo "         ${GENOMIC_FILE}"
  echo "         Run run_cohort_analysis_pipeline.sh (step 1) first."
fi

if [[ ! -f "${CLINICAL_FILE}" ]]; then
  echo "  [WARN] Clinical features file not found:"
  echo "         ${CLINICAL_FILE}"
  echo "         Run run_cohort_analysis_pipeline.sh (step 2) first."
fi

if [[ ! -f "${RAW_FEATURES_FILE}" ]]; then
  echo "  ERROR: Raw features file not found: ${RAW_FEATURES_FILE}"
  exit 1
fi

PRS_FILE_COUNT=$(find "${SCORES}/prsScores" -name "*mixed*.csv" 2>/dev/null | wc -l | tr -d ' ')
if [[ "${PRS_FILE_COUNT}" -eq 0 ]]; then
  echo "  ERROR: No *.mixed.prs.csv files found in ${SCORES}/prsScores/"
  exit 1
fi
echo "  Found ${PRS_FILE_COUNT} *.mixed.prs.csv files"

# Check at least one of genomic/clinical prerequisite exists
if [[ ! -f "${GENOMIC_FILE}" && ! -f "${CLINICAL_FILE}" ]]; then
  echo "  ERROR: Neither genomic nor clinical cohort analysis outputs found."
  echo "         Please run run_cohort_analysis_pipeline.sh first."
  exit 1
fi

echo ""
echo "Running cohort bNMF..."
echo ""

# ── Build optional argument list ──────────────────────────────────────────────
EXTRA_ARGS=()

if [[ -n "${COHORTS:-}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS+=(--cohorts ${COHORTS})
fi

if [[ "${NO_HOLDOUT_PRS}" == "1" ]]; then
  EXTRA_ARGS+=(--no_holdout_prs)
fi

# ── Run script ────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")/../scripts" && pwd)"

python "${SCRIPT_DIR}/cohort_bnmf_subtype_analysis.py" \
  --pheno_data             "${PHENO_DATA}" \
  --raw_features_file      "${RAW_FEATURES_FILE}" \
  --mode                   "${MODE}" \
  --filter_strategy        "${FILTER_STRATEGY}" \
  --use_set                "${USE_SET}" \
  --k_min                  "${K_MIN}" \
  --k_max                  "${K_MAX}" \
  --n_runs                 "${N_RUNS}" \
  --k_select_method        "${K_SELECT_METHOD}" \
  --l1_ratio               "${L1_RATIO}" \
  --alpha_W                "${ALPHA_W}" \
  --min_effect_size        "${MIN_EFFECT_SIZE}" \
  --specificity_tier_max   "${SPECIFICITY_TIER_MAX}" \
  --high_risk_threshold    "${HIGH_RISK_THRESHOLD}" \
  --low_control_threshold  "${LOW_CONTROL_THRESHOLD}" \
  "${EXTRA_ARGS[@]}"

echo ""
echo "============================================================"
echo "Outputs"
echo "============================================================"
if [[ "${MODE}" == "per_cohort" || "${MODE}" == "both" ]]; then
  echo "  Per-cohort results   : ${SCORES}/cohortBnmf/{cohort}/"
  echo "    cluster_assignments.csv   — soft + hard cluster membership per IID"
  echo "    feature_loadings.csv      — H matrix (feature weights per cluster)"
  echo "    basis_matrix.csv          — normalised W membership columns"
  echo "    cophenetic_curve.csv      — cophenetic correlation per k"
  echo "    cluster_profile.csv       — per-cluster median statistics"
  echo "  Overview             : ${SCORES}/cohortBnmf/cohort_bnmf_overview.csv"
  echo "  Figures              : ${PHENO_DATA}/figures/cohortBnmf/{cohort}/"
  echo ""
fi
if [[ "${MODE}" == "combined" || "${MODE}" == "both" ]]; then
  echo "  Combined results     : ${SCORES}/combinedCohortBnmf/"
  echo "    cluster_assignments.csv      — membership + is_high_risk_{cohort} flags"
  echo "    feature_loadings.csv         — H matrix ({cohort}_clin/gen_ prefixed columns)"
  echo "    cohort_cluster_affinity.csv  — cohort × cluster proportions of high-risk individuals"
  echo "    weighted_feature_list.csv    — per-cluster features ranked by loading × cohort weight"
  echo "    cluster_dominant_cohort.csv  — which cohort dominates each cluster"
  echo "  Figures              : ${PHENO_DATA}/figures/combinedCohortBnmf/"
  echo "    cohort_cluster_affinity.png  — heatmap of cohort–cluster affinities"
  echo ""
fi
echo "Done."
