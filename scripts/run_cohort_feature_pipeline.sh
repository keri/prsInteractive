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
#   Step 5  permutation_test_cluster_enrichment.py
#   Step 6  cluster_subtype_dotmatrix.py          (Figure 8)
#   Step 7  subtype_cohort_stackedbars.py         (companion figure)
#
#   The circos and donut-pie steps are retired; see the block at the end.
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
#   --run-exclusive             Run ExclusiveHighCases bootstrap in Step 1 (slow; adds --run_exclusive flag)
#   --exclusive-iterations INT  Bootstrap iterations for ExclusiveHighCases (default: 1000)
#   --skip-calculate            Skip Step 1 (use existing importantFeaturesAcrossCohortsAndTrainingData.Filtered.*.csv)
#   --n-perm             INT    Permutations for Step 5 (default: 10000)
#   --perm-seed          INT    Seed for Step 5 (default: 42)
#   --only-plots                Skip Steps 1-4, regenerate the test and figures only
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
N_PERM="${N_PERM:-10000}"
PERM_SEED="${PERM_SEED:-42}"

SKIP_CALCULATE=0
ONLY_PLOTS=0
RUN_EXCLUSIVE=0
EXCLUSIVE_ITERATIONS=1000

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
        --n-perm)            N_PERM="$2";              shift 2 ;;
        --perm-seed)         PERM_SEED="$2";           shift 2 ;;
        --run-exclusive)     RUN_EXCLUSIVE=1;            shift ;;
        --exclusive-iterations) EXCLUSIVE_ITERATIONS="$2"; shift 2 ;;
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
PERM_CLUSTER_CSV="${ENRICHMENT_OUT}/permutation_enrichment_cluster_key.csv"

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
echo "  run_exclusive    : ${RUN_EXCLUSIVE}"
echo "  excl_iterations  : ${EXCLUSIVE_ITERATIONS}"
echo "  n_perm           : ${N_PERM}"
echo "  perm_seed        : ${PERM_SEED}"
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

    CALC_ARGS=(
        --pheno_data              "${PHENO_DATA}"
        --feature_scores_file_filtered "${FEATURE_SCORES}"
        --filter_strategy         "${FILTER_STRATEGY}"
    )
    if [[ ${RUN_EXCLUSIVE} -eq 1 ]]; then
        CALC_ARGS+=( --run_exclusive --n_exclusive_iterations "${EXCLUSIVE_ITERATIONS}" )
    fi

    python "${SCRIPT_DIR}/calculate_top_features_in_cohort.py" "${CALC_ARGS[@]}"

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
#  STEP 5 — Permutation test for differential cluster/subtype correspondence
# =============================================================================
step_header "Permutation enrichment test (permutation_test_cluster_enrichment.py)"

python "${SCRIPT_DIR}/permutation_test_cluster_enrichment.py" \
    --assignments_file "${ASSIGNMENTS_CSV}" \
    --out_dir          "${ENRICHMENT_OUT}" \
    --n_perm           "${N_PERM}" \
    --seed             "${PERM_SEED}"

echo "  Output: ${PERM_CLUSTER_CSV}"

# =============================================================================
#  STEP 6 — Figure 8: cluster x subtype dot matrix
# =============================================================================
step_header "Cluster-subtype dot matrix (cluster_subtype_dotmatrix.py)"

mkdir -p "${PLOTS_OUT}"

python "${SCRIPT_DIR}/cluster_subtype_dotmatrix.py" \
    --data            "${ASSIGNMENTS_CSV}" \
    --permutation_csv "${PERM_CLUSTER_CSV}" \
    --out_dir         "${PLOTS_OUT}"

echo "  Output: ${PLOTS_OUT}/cluster_subtype_dotmatrix.{png,pdf}"

# =============================================================================
#  STEP 7 — Companion stacked-bar figure
# =============================================================================
step_header "Subtype cohort stacked bars (subtype_cohort_stackedbars.py)"

python "${SCRIPT_DIR}/subtype_cohort_stackedbars.py" \
    --assignments "${ASSIGNMENTS_CSV}" \
    --out_dir     "${PLOTS_OUT}"

echo "  Output: ${PLOTS_OUT}/subtype_cohort_bars_*.png"

# =============================================================================
#  RETIRED PLOT STEPS — circos + donut pies
# =============================================================================
# Neither figure is used in the manuscript.  Both also predate the move to
# signed concordance and would need reworking before they could be trusted:
# circos_subtype_cohort.py reads the sign of effect_size_r directly (bar
# outward/inward) and subtype_cohort_pies.py applies the old
# {+1: 1.0, 0: 0.5, -1: 0.0} weighting, which drops discordant features.
# Re-enable only after updating them for the signed convention.
#
# step_header "Circos plots (circos_subtype_cohort.py)"
# python "${SCRIPT_DIR}/circos_subtype_cohort.py" \
#     --assignments "${ASSIGNMENTS_CSV}" \
#     --out_dir     "${PLOTS_OUT}"
#
# step_header "Subtype cohort pie charts (subtype_cohort_pies.py)"
# python "${SCRIPT_DIR}/subtype_cohort_pies.py" \
#     --assignments "${ASSIGNMENTS_CSV}" \
#     --out_dir     "${PLOTS_OUT}"

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Pipeline complete."
echo "  Enrichment results : ${ENRICHMENT_OUT}"
echo "  Plots              : ${PLOTS_OUT}"
echo "════════════════════════════════════════════════════════════════"
