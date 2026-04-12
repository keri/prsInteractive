#!/usr/bin/env bash
# =============================================================================
# run_cohort_feature_pipeline.sh
#
# End-to-end pipeline for the PRS cohort important-feature analysis:
#
#   Step 1  calculate_top_features_in_cohort.py   [skippable with --skip-calculate]
#   Step 2  enrich_exclusive_cohort_features.py
#   Step 3  compare_clinical_features_across_cohorts.py
#   Step 4  clusterSubtypeAssignment.py
#   Step 5  circos_subtype_cohort.py
#   Step 6  subtype_cohort_pies.py
#
# Usage
# -----
#   bash scripts/run_cohort_feature_pipeline.sh [OPTIONS]
#
# Options
#   --pheno-data         PATH   Path to combinedAnalysis directory
#                               (default: results/type2Diabetes/combinedAnalysis)
#   --subtype-dir        PATH   Path to subtypeResearch directory
#                               (default: results/type2Diabetes/subtypeResearch)
#   --feature-scores     PATH   featureScoresReducedFinalModel.combined.filtered.csv
#                               (default: <pheno-data>/scores/featureScoresReducedFinalModel.combined.filtered.csv)
#   --raw-features       PATH   Raw participant environment CSV for clinical MWU
#                               (default: results/participant_environment.csv)
#   --use-set            STR    holdout | validation | all  (default: holdout)
#   --filter-strategy    STR    tiered | strict | majority | differential  (default: tiered)
#   --balance-strategy   STR    unbalanced | balanced  (default: unbalanced)
#   --skip-calculate            Skip Step 1 (use existing importantFeaturesAcrossCohortsAndTrainingData.Filtered.*.csv)
#   --only-plots                Skip Steps 1-4, regenerate circos + pie charts only
#   --help
#
# Environment variables (override defaults)
#   PHENO_DATA, SUBTYPE_DIR, FEATURE_SCORES_FILE_FILTERED, RAW_FEATURE_SCORES,
#   USE_SET, FILTER_STRATEGY, BALANCE_STRATEGY
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Defaults ─────────────────────────────────────────────────────────────────
PHENO_DATA="${PHENO_DATA:-${REPO_ROOT}/results/type2Diabetes/combinedAnalysis}"
SUBTYPE_DIR="${SUBTYPE_DIR:-${REPO_ROOT}/results/type2Diabetes/subtypeResearch}"
RAW_FEATURE_SCORES="${RAW_FEATURE_SCORES:-${REPO_ROOT}/results/participant_environment.csv}"
USE_SET="${USE_SET:-holdout}"
FILTER_STRATEGY="${FILTER_STRATEGY:-tiered}"
BALANCE_STRATEGY="${BALANCE_STRATEGY:-unbalanced}"

SKIP_CALCULATE=0
ONLY_PLOTS=0

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pheno-data)        PHENO_DATA="$2";          shift 2 ;;
        --subtype-dir)       SUBTYPE_DIR="$2";         shift 2 ;;
        --feature-scores)    FEATURE_SCORES="$2";      shift 2 ;;
        --raw-features)      RAW_FEATURE_SCORES="$2";  shift 2 ;;
        --use-set)           USE_SET="$2";              shift 2 ;;
        --filter-strategy)   FILTER_STRATEGY="$2";     shift 2 ;;
        --balance-strategy)  BALANCE_STRATEGY="$2";    shift 2 ;;
        --skip-calculate)    SKIP_CALCULATE=1;          shift ;;
        --only-plots)        ONLY_PLOTS=1;              shift ;;
        --help|-h)
            sed -n '/^# Usage/,/^# ====/p' "$0" | grep -v '^# ====' | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Default feature scores path if not explicitly set
FEATURE_SCORES="${FEATURE_SCORES:-${PHENO_DATA}/scores/featureScoresReducedFinalModel.combined.filtered.csv}"

# Derived paths
SCORES_PATH="${PHENO_DATA}/scores"
TIERED_CSV="${SCORES_PATH}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.Filtered.${FILTER_STRATEGY}.csv"
ENRICHMENT_CSV="${SCORES_PATH}/functionalEnrichment/functional_cluster_enrichment_${FILTER_STRATEGY}.csv"
ASSIGNMENTS_CSV="${SUBTYPE_DIR}/enrichment/subtype_feature_assignments.csv"
CLUSTER_SUBTYPE_MAP="${SUBTYPE_DIR}/cluster_subtype_map.csv"
ENRICHMENT_OUT="${SUBTYPE_DIR}/enrichment"
PLOTS_OUT="${ENRICHMENT_OUT}/plots"

# ── Helpers ───────────────────────────────────────────────────────────────────
STEP=0
step_header() {
    STEP=$((STEP + 1))
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Step ${STEP}: $1"
    echo "════════════════════════════════════════════════════════════════"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: Required file not found: $1" >&2
        exit 1
    fi
}

# ── Pre-flight checks ─────────────────────────────────────────────────────────
echo "PRS Cohort Feature Pipeline"
echo "==========================="
echo "  pheno_data       : ${PHENO_DATA}"
echo "  subtype_dir      : ${SUBTYPE_DIR}"
echo "  feature_scores   : ${FEATURE_SCORES}"
echo "  raw_features     : ${RAW_FEATURE_SCORES}"
echo "  use_set          : ${USE_SET}"
echo "  filter_strategy  : ${FILTER_STRATEGY}"
echo "  balance_strategy : ${BALANCE_STRATEGY}"
echo "  skip_calculate   : ${SKIP_CALCULATE}"
echo "  only_plots       : ${ONLY_PLOTS}"
echo ""

if [[ ${ONLY_PLOTS} -eq 0 ]]; then
    check_file "${FEATURE_SCORES}"
    check_file "${RAW_FEATURE_SCORES}"
    check_file "${CLUSTER_SUBTYPE_MAP}"
fi

# =============================================================================
#  STEP 1 — SHAP feature importance + cross-cohort filtering
# =============================================================================
if [[ ${ONLY_PLOTS} -eq 0 && ${SKIP_CALCULATE} -eq 0 ]]; then
    step_header "SHAP feature importance + cross-cohort filtering (calculate_top_features_in_cohort.py)"

    python "${SCRIPT_DIR}/calculate_top_features_in_cohort.py" \
        --pheno_data              "${PHENO_DATA}" \
        --feature_scores_file_filtered "${FEATURE_SCORES}" \
        --filter_strategy         "${FILTER_STRATEGY}"

    check_file "${TIERED_CSV}"
    echo "  Output: ${TIERED_CSV}"

elif [[ ${ONLY_PLOTS} -eq 0 && ${SKIP_CALCULATE} -eq 1 ]]; then
    step_header "SHAP feature importance — SKIPPED (--skip-calculate)"
    check_file "${TIERED_CSV}"
    echo "  Using existing: ${TIERED_CSV}"
fi

# =============================================================================
#  STEP 2 — Functional cluster enrichment + clinical MWU
# =============================================================================
if [[ ${ONLY_PLOTS} -eq 0 ]]; then
    step_header "Functional cluster enrichment + clinical MWU (enrich_exclusive_cohort_features.py)"

    python "${SCRIPT_DIR}/enrich_exclusive_cohort_features.py" \
        --pheno_data          "${PHENO_DATA}" \
        --raw_feature_scores  "${RAW_FEATURE_SCORES}" \
        --filter_strategy     "${FILTER_STRATEGY}" \
        --use_set             "${USE_SET}" \
        --balance_strategy    "${BALANCE_STRATEGY}"

    check_file "${ENRICHMENT_CSV}"
    echo "  Output: ${ENRICHMENT_CSV}"
fi

# =============================================================================
#  STEP 3 — Cross-cohort clinical comparison (all three risk tiers)
# =============================================================================
if [[ ${ONLY_PLOTS} -eq 0 ]]; then
    step_header "Cross-cohort clinical comparison (compare_clinical_features_across_cohorts.py)"

    python "${SCRIPT_DIR}/compare_clinical_features_across_cohorts.py" \
        --pheno_data         "${PHENO_DATA}" \
        --raw_features_file  "${RAW_FEATURE_SCORES}" \
        --use_set            "${USE_SET}"

    echo "  Output: ${SCORES_PATH}/cohortClinicalComparison/"
fi

# =============================================================================
#  STEP 4 — Cluster-subtype assignment
# =============================================================================
if [[ ${ONLY_PLOTS} -eq 0 ]]; then
    step_header "Cluster-subtype assignment (clusterSubtypeAssignment.py)"

    mkdir -p "${ENRICHMENT_OUT}"

    python "${SCRIPT_DIR}/clusterSubtypeAssignment.py" \
        --tiered_csv          "${ENRICHMENT_CSV}" \
        --cluster_subtype_map "${CLUSTER_SUBTYPE_MAP}" \
        --out_dir             "${ENRICHMENT_OUT}"

    check_file "${ASSIGNMENTS_CSV}"
    echo "  Output: ${ASSIGNMENTS_CSV}"
fi

# =============================================================================
#  STEP 5 — Circos plots (4 variants per cohort)
# =============================================================================
step_header "Circos plots (circos_subtype_cohort.py)"

mkdir -p "${PLOTS_OUT}"

python "${SCRIPT_DIR}/circos_subtype_cohort.py" \
    --assignments "${ASSIGNMENTS_CSV}" \
    --out_dir     "${PLOTS_OUT}"

echo "  Output: ${PLOTS_OUT}/circos_*.png"

# =============================================================================
#  STEP 6 — Subtype donut pie charts
# =============================================================================
step_header "Subtype cohort pie charts (subtype_cohort_pies.py)"

python "${SCRIPT_DIR}/subtype_cohort_pies.py" \
    --assignments "${ASSIGNMENTS_CSV}" \
    --out_dir     "${PLOTS_OUT}"

echo "  Output: ${PLOTS_OUT}/subtype_cohort_pies_*.png"

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Pipeline complete."
echo "  Enrichment results : ${ENRICHMENT_OUT}"
echo "  Plots              : ${PLOTS_OUT}"
echo "════════════════════════════════════════════════════════════════"
