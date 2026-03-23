#!/usr/bin/env python3

#!/usr/bin/env bash
# run_gxg_gplus_analysis.sh
#
# Runs combine_gxg_gplus_prs_analysis.py (Step 1–7 pipeline), then on
# success runs plot_prs_composite_comparison_to_main.py.
#
# Usage:
#   bash run_gxg_gplus_analysis.sh [options]
#
# Required:
#   -p  PHENO            Phenotype label (e.g. type2Diabetes)
#   -r  PHENO_RESULTS    Root results dir containing productEpi/ and summedEpi/
#   -R  RESULTS_PATH     Broader results dir (used to locate participant_environment.csv)
#
# Optional:
#   -f  FEATURE_SCORES   Path to pre-built combined feature-scores CSV
#   -e  RAW_ENV_FILE     Path to raw environmental measures CSV
#   -s  STEPS            Steps to run (default: 1,2,3,4,5,6,7)
#   -t  THRESHOLD        SHAP z-score threshold (default: 1.99)
#   -a                   Include prs_all models (--use_all flag)
#   -S  SCRIPT_DIR       Directory containing the Python scripts
#                        (default: same directory as this shell script)

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
PHENO=""
PHENO_RESULTS=""
RESULTS_PATH=""
FEATURE_SCORES=""
RAW_ENV_FILE=""
STEPS="1,2,3,4,5,6,7"
THRESHOLD="1.99"
USE_ALL=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COMBINE_SCRIPT="combine_gxg_gplus_prs_analysis.py"
PLOT_SCRIPT="plot_prs_composite_comparison_to_main.py"

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,2\}//'
    exit 1
}

while getopts ":p:r:R:f:e:s:t:S:ah" opt; do
    case $opt in
        p) PHENO="$OPTARG" ;;
        r) PHENO_RESULTS="$OPTARG" ;;
        R) RESULTS_PATH="$OPTARG" ;;
        f) FEATURE_SCORES="$OPTARG" ;;
        e) RAW_ENV_FILE="$OPTARG" ;;
        s) STEPS="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        a) USE_ALL="--use_all" ;;
        S) SCRIPT_DIR="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: -$OPTARG requires an argument." >&2; exit 1 ;;
        \?) echo "ERROR: Unknown option -$OPTARG." >&2; exit 1 ;;
    esac
done

# ── Validate required args ────────────────────────────────────────────────────
if [[ -z "$PHENO" || -z "$PHENO_RESULTS" || -z "$RESULTS_PATH" ]]; then
    echo "ERROR: -p (PHENO), -r (PHENO_RESULTS), and -R (RESULTS_PATH) are required." >&2
    echo "Run with -h for usage." >&2
    exit 1
fi

# ── Derived paths ─────────────────────────────────────────────────────────────
# The combined analysis output directory that Step 1 produces.
COMBINED_ANALYSIS_DIR="${PHENO_RESULTS}/combinedAnalysis"

# ── Optional argument assembly ────────────────────────────────────────────────
FEATURE_SCORES_ARG=""
[[ -n "$FEATURE_SCORES" ]] && FEATURE_SCORES_ARG="--feature_scores_file ${FEATURE_SCORES}"

RAW_ENV_ARG=""
[[ -n "$RAW_ENV_FILE" ]] && RAW_ENV_ARG="--raw_env_features_file ${RAW_ENV_FILE}"

# ── Step 1: combine_gxg_gplus_prs_analysis.py ────────────────────────────────
echo ""
echo "============================================================"
echo "  combine_gxg_gplus_prs_analysis.py"
echo "  phenotype    : ${PHENO}"
echo "  pheno_results: ${PHENO_RESULTS}"
echo "  results_path : ${RESULTS_PATH}"
echo "  steps        : ${STEPS}"
echo "============================================================"
echo ""

python "${SCRIPT_DIR}/${COMBINE_SCRIPT}" \
    --pheno_results "${PHENO_RESULTS}" \
    --pheno         "${PHENO}" \
    --results_path  "${RESULTS_PATH}" \
    --steps         "${STEPS}" \
    --threshold     "${THRESHOLD}" \
    ${FEATURE_SCORES_ARG} \
    ${RAW_ENV_ARG} \
    ${USE_ALL}
    
echo ""
echo "  combine_gxg_gplus_prs_analysis.py completed successfully."
echo ""

# ── Step 2: plot_prs_composite_comparison_to_main.py ─────────────────────────
echo "============================================================"
echo "  plot_prs_composite_comparison_to_main.py"
echo "  PHENO_DATA: ${COMBINED_ANALYSIS_DIR}"
echo "============================================================"
echo ""

export PHENO_DATA="${COMBINED_ANALYSIS_DIR}"

python "${SCRIPT_DIR}/${PLOT_SCRIPT}"

echo ""
echo "  plot_prs_composite_comparison_to_main.py completed successfully."
echo ""
echo "============================================================"
echo "  All done."
echo "  Combined analysis : ${COMBINED_ANALYSIS_DIR}"
echo "============================================================"