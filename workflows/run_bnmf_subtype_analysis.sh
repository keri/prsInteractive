#!/bin/bash
# =============================================================================
# run_bnmf_subtype_analysis.sh
#
# Workflow: bNMF phenotype subtype discovery
#
# Runs prs_bnmf_subtype_analysis.py to identify phenotypic subtypes by
# applying Bayesian-inspired Non-negative Matrix Factorization (bNMF) to
# ePRS model scores combined with clinical measures.
#
# Inspired by Udler et al. 2018 PLOS Medicine T2D subtype discovery.
# DOI: https://doi.org/10.1371/journal.pmed.1002654
#
# Usage:
#   bash run_bnmf_subtype_analysis.sh <PHENO> [EPI_COMBO] [MODE] [OPTIONS]
#
# Arguments:
#   PHENO      Phenotype name (e.g. type2Diabetes, celiacDisease)
#   EPI_COMBO  Epistatic combination: "sum" (default) or "prod"
#   MODE       Input mode: "prs_scores" (default) or "weighted_features"
#
# Optional environment overrides (export before running):
#   BNMF_K_MIN              Minimum k to test (default: 2)
#   BNMF_K_MAX              Maximum k to test (default: 8)
#   BNMF_N_RUNS             NMF restarts per k  (default: 30)
#   BNMF_CASES_ONLY         Set to "true" to restrict to cases only
#   BNMF_K_SELECT_METHOD    "elbow" or "max"  (default: elbow)
#   BNMF_L1_RATIO           L1 regularization ratio  (default: 0.1)
#   BNMF_ALPHA_W            Regularization strength on W  (default: 0.1)
#
# Example (T2D, summed epistasis, PRS-score mode):
#   bash run_bnmf_subtype_analysis.sh type2Diabetes sum prs_scores
#
# Example (T2D, weighted features, cases only):
#   BNMF_CASES_ONLY=true \
#   bash run_bnmf_subtype_analysis.sh type2Diabetes sum weighted_features
# =============================================================================

set -euo pipefail

PHENO="${1:?Usage: $0 <PHENO> [EPI_COMBO] [MODE]}"
EPI_COMBO="${2:-sum}"
MODE="${3:-prs_scores}"

# ------------------------------------------------------------------
# Source global environment config
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_CONFIG="${SCRIPT_DIR}/../env.config"

if [ ! -f "${ENV_CONFIG}" ]; then
    echo "ERROR: env.config not found at ${ENV_CONFIG}"
    exit 1
fi
source "${ENV_CONFIG}"

# ------------------------------------------------------------------
# Resolve phenotype-specific paths
# ------------------------------------------------------------------
if [ "${EPI_COMBO}" == "sum" ]; then
    COMBO_FOLDER='summedEpi'
else
    COMBO_FOLDER='productEpi'
fi

PHENO_DATA="${RESULTS_PATH}/${PHENO}/${COMBO_FOLDER}"
CONFIG_FILE="${PHENO_DATA}/pheno.config"

if [ ! -d "${PHENO_DATA}" ]; then
    echo "ERROR: Phenotype directory not found: ${PHENO_DATA}"
    echo "Run envSetUp.sh first to create the phenotype environment."
    exit 1
fi

source "${CONFIG_FILE}"

# ------------------------------------------------------------------
# Resolve clinical measures file
# ------------------------------------------------------------------
CLINICAL_FILE="${ENV_FILE:-}"
if [ -z "${CLINICAL_FILE}" ] || [ ! -f "${CLINICAL_FILE}" ]; then
    # Fallback: look for participant_environment.csv in data directory
    CLINICAL_FILE="${DATA_PATH}/participant_environment.csv"
fi
if [ ! -f "${CLINICAL_FILE}" ]; then
    echo "ERROR: Clinical measures file not found."
    echo "Set ENV_FILE in pheno.config or ensure ${DATA_PATH}/participant_environment.csv exists."
    exit 1
fi

# ------------------------------------------------------------------
# bNMF parameters (with environment variable overrides)
# ------------------------------------------------------------------
K_MIN="${BNMF_K_MIN:-2}"
K_MAX="${BNMF_K_MAX:-8}"
N_RUNS="${BNMF_N_RUNS:-30}"
K_SELECT_METHOD="${BNMF_K_SELECT_METHOD:-elbow}"
L1_RATIO="${BNMF_L1_RATIO:-0.1}"
ALPHA_W="${BNMF_ALPHA_W:-0.1}"
CASES_ONLY_FLAG=""
if [ "${BNMF_CASES_ONLY:-false}" == "true" ]; then
    CASES_ONLY_FLAG="--cases_only"
fi

# ------------------------------------------------------------------
# Important features file (optional — for weighted_features mode)
# ------------------------------------------------------------------
IMPORTANT_FEATURES_FILE="${PHENO_DATA}/scores/featureScoresReducedFinalModel.filtered.csv"
IMPORTANT_FEATURES_ARG=""
if [ -f "${IMPORTANT_FEATURES_FILE}" ]; then
    IMPORTANT_FEATURES_ARG="--important_features_file ${IMPORTANT_FEATURES_FILE}"
    echo "  Using important features from: ${IMPORTANT_FEATURES_FILE}"
else
    echo "  No filtered feature scores found at ${IMPORTANT_FEATURES_FILE}."
    echo "  In weighted_features mode all features will be used (may be slow)."
fi

# ------------------------------------------------------------------
# Validate prerequisites
# ------------------------------------------------------------------
echo ""
echo "======================================================="
echo "  bNMF Subtype Analysis"
echo "  Phenotype : ${PHENO}"
echo "  EpiCombo  : ${COMBO_FOLDER}"
echo "  Mode      : ${MODE}"
echo "  k range   : ${K_MIN}–${K_MAX}  (${N_RUNS} runs each)"
echo "  k select  : ${K_SELECT_METHOD}"
echo "  PhenoData : ${PHENO_DATA}"
echo "  Clinical  : ${CLINICAL_FILE}"
echo "======================================================="
echo ""

REQUIRED_SCRIPT="${SCRIPTS_DIR}/prs_bnmf_subtype_analysis.py"
if [ ! -f "${REQUIRED_SCRIPT}" ]; then
    echo "ERROR: Script not found: ${REQUIRED_SCRIPT}"
    exit 1
fi

# Check that at least one PRS file exists
PRS_SCORES_DIR="${PHENO_DATA}/scores/prsScores"
if [ ! -d "${PRS_SCORES_DIR}" ]; then
    echo "ERROR: PRS scores directory not found: ${PRS_SCORES_DIR}"
    echo "Run run_prs_calculations.sh first."
    exit 1
fi
PRS_COUNT=$(find "${PRS_SCORES_DIR}" -name "*mixed*.csv" ! -name "*holdout*" | wc -l | tr -d ' ')
if [ "${PRS_COUNT}" -eq 0 ]; then
    echo "ERROR: No *.mixed.prs.csv files found in ${PRS_SCORES_DIR}."
    echo "Run run_prs_calculations.sh first."
    exit 1
fi
echo "  Found ${PRS_COUNT} PRS model file(s) in ${PRS_SCORES_DIR}"

# ------------------------------------------------------------------
# Export variables for Python script
# ------------------------------------------------------------------
export PHENO_DATA
export PHENO

# ------------------------------------------------------------------
# Run bNMF analysis
# ------------------------------------------------------------------
echo ""
echo "[RUNNING] prs_bnmf_subtype_analysis.py ..."
echo ""

python "${REQUIRED_SCRIPT}" \
    --pheno_data "${PHENO_DATA}" \
    --clinical_file "${CLINICAL_FILE}" \
    --mode "${MODE}" \
    --k_min "${K_MIN}" \
    --k_max "${K_MAX}" \
    --n_runs "${N_RUNS}" \
    --l1_ratio "${L1_RATIO}" \
    --alpha_w "${ALPHA_W}" \
    --k_select_method "${K_SELECT_METHOD}" \
    --phenotype_name "${PHENO}" \
    ${IMPORTANT_FEATURES_ARG} \
    ${CASES_ONLY_FLAG}

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo ""
    echo "ERROR: prs_bnmf_subtype_analysis.py failed (exit code ${EXIT_CODE})"
    exit ${EXIT_CODE}
fi

# ------------------------------------------------------------------
# Report outputs
# ------------------------------------------------------------------
OUT_DIR="${PHENO_DATA}/scores/bnmf"
FIG_DIR="${PHENO_DATA}/figures/bnmf"

echo ""
echo "======================================================="
echo "  COMPLETED SUCCESSFULLY"
echo ""
echo "  Output files:"
for f in \
    cluster_assignments.csv \
    feature_loadings.csv \
    cophenetic_curve.csv \
    cluster_clinical_profiles.csv \
    prs_clinical_associations.csv \
    udler_similarity.csv; do
    TARGET="${OUT_DIR}/${f}"
    if [ -f "${TARGET}" ]; then
        echo "    ✓ ${TARGET}"
    fi
done

echo ""
echo "  Figures:"
for f in "${FIG_DIR}"/*.png; do
    [ -f "${f}" ] && echo "    ✓ ${f}"
done
echo "======================================================="
