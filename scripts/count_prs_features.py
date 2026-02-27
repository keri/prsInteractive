"""
count_prs_features.py
---------------------
Count unique G, GxG, and E features used in cardio_product and cardio_summed
PRS calculations from featureScoresReducedFinalModel.combined.filtered.csv.

Only rows with:
  - model  != 'covar'   (non-covariate)
  - coefs  != 0         (non-zero beta coefficient)
are counted.

Feature structure (positional, comma-separated)
------------------------------------------------
The 'feature' column encodes up to three components per cell.
Classification is purely positional — no prefix matching needed:

  Token count  │ First token  │ Remaining  │ Components extracted
  ─────────────┼──────────────┼────────────┼──────────────────────
  1            │ rsID         │ —          │ G   = token[0]
  2            │ rsID         │ rsID       │ GxG = (token[0], token[1])
  2            │ non-rsID     │ rsID       │ E   = token[0]
               │              │            │ G   = token[1]
  3            │ non-rsID     │ rsID, rsID │ E   = token[0]
               │              │            │ GxG = (token[1], token[2])

A single feature cell can therefore contribute to more than one type counter
(e.g. an E×GxG row increments both the E set and the GxG set).

Output files written to {pheno_data}/scores/
  feature_counts_summary.csv   – wide table: one row per model
  feature_counts_long.csv      – long table: one row per model × feature_type
  feature_lists/
      {model}_G_features.csv   – unique G rsIDs per model
      {model}_GxG_features.csv – unique GxG pairs per model  (as "rs1,rs2")
      {model}_E_features.csv   – unique E features per model

Usage
-----
python count_prs_features.py \\
    --pheno_data  /path/to/pheno_data \\
    --input_file  featureScoresReducedFinalModel.combined.filtered.csv
"""

import logging
import argparse
from pathlib import Path

import pandas as pd
import warnings
warnings.filterwarnings("ignore")

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
PRS_MODELS  = ("cardio_product", "cardio_summed")
COEF_COL    = "coefs"
MODEL_COL   = "model"
FEATURE_COL = "feature"
COVAR_LABEL = "covariate_summed"


# ---------------------------------------------------------------------------
# Feature parsing
# ---------------------------------------------------------------------------

def _is_rsid(token: str) -> bool:
    """True if *token* looks like an rsID (starts with 'rs' followed by digits)."""
    t = token.strip().lower()
    return t.startswith("rs") and t[2:].isdigit()


def _parse_feature_components(feature: str) -> dict[str, str | None]:
    """
    Parse one 'feature' cell and return its component parts.

    Returns
    -------
    {
        'E'  : str | None   – E feature token (non-rsID first token)
        'G'  : str | None   – single rsID when feature is G or E×G
        'GxG': str | None   – canonical "rs1,rs2" string when GxG or E×GxG
    }

    All values are lower-cased for consistent deduplication.
    """
    tokens = [t.strip().lower() for t in str(feature).split(",") if t.strip()]
    out    = {"E": None, "G": None, "GxG": None}

    if not tokens:
        return out

    first_is_rsid = _is_rsid(tokens[0])

    # ── Pure G : single rsID ────────────────────────────────────────────────
    if len(tokens) == 1:
        if first_is_rsid:
            out["G"] = tokens[0]
        # else: single non-rsID token → covariate/other, all None
        return out

    # ── Two tokens ──────────────────────────────────────────────────────────
    if len(tokens) == 2:
        if first_is_rsid and _is_rsid(tokens[1]):
            # Pure GxG: rs1,rs2
            out["GxG"] = f"{tokens[0]},{tokens[1]}"
        elif not first_is_rsid and _is_rsid(tokens[1]):
            # E × G: e_feature,rs1
            out["E"] = tokens[0]
            out["G"] = tokens[1]
        else:
            log.debug("Unrecognised 2-token feature: %s", feature)
        return out

    # ── Three tokens ────────────────────────────────────────────────────────
    if len(tokens) == 3:
        if not first_is_rsid and _is_rsid(tokens[1]) and _is_rsid(tokens[2]):
            # E × GxG: e_feature,rs1,rs2
            out["E"]   = tokens[0]
            out["GxG"] = f"{tokens[1]},{tokens[2]}"
        else:
            log.debug("Unrecognised 3-token feature: %s", feature)
        return out

    # ── More than 3 tokens — unexpected ─────────────────────────────────────
    log.warning("Feature with >3 tokens (skipped): %s", feature)
    return out


# ---------------------------------------------------------------------------
# Core counting logic
# ---------------------------------------------------------------------------

def count_features(df: pd.DataFrame) -> dict[str, dict[str, object]]:
    """
    For each PRS model in PRS_MODELS, build deduplicated sets of
    G, GxG, and E components from non-covariate, non-zero-coef rows.

    Covariate exclusion strategy
    ----------------------------
    Features that appear in ANY row where model == 'covariate_summed' are collected into
    a normalised exclusion set.  That set is then used to drop matching
    features from the target model's rows, regardless of what model label
    those rows carry.  Both sides are stripped and lower-cased before
    comparison so whitespace/case differences never cause a miss.

    Returns
    -------
    {
      'cardio_product': {
          'G'  : {'count': 123, 'features': {'rs123', ...}},
          'GxG': {'count':  45, 'features': {'rs1,rs2', ...}},
          'E'  : {'count':  67, 'features': {'age', ...}},
      },
      ...
    }
    """
    # Build normalised covariate feature exclusion set once, used for all models
    covar_features: set[str] = set(
        df.loc[
            df[MODEL_COL].str.strip().str.lower() == COVAR_LABEL,
            FEATURE_COL,
        ]
        .dropna()
        .str.strip()
        .str.lower()
    )
    log.info("Covariate features to exclude: %d", len(covar_features))

    # Normalised feature column used only for the exclusion mask
    feature_normalised = df[FEATURE_COL].str.strip().str.lower()

    results: dict[str, dict] = {}

    for model in PRS_MODELS:
        log.info("Processing model: %s", model)

        mask = (
            (df[MODEL_COL].str.strip().str.lower() == model.lower()) &
            (~feature_normalised.isin(covar_features))                &
            (df[COEF_COL].fillna(0) != 0)
        )
        subset = df[mask].copy()
        log.info("  Rows after filtering: %d", len(subset))

        if subset.empty:
            log.warning("  No rows found for '%s' after filtering.", model)
            results[model] = {ft: {"count": 0, "features": set()} for ft in ("G", "GxG", "E")}
            continue

        # Accumulate components across all rows
        g_set:   set[str] = set()
        gxg_set: set[str] = set()
        e_set:   set[str] = set()

        for feat in subset[FEATURE_COL].dropna():
            components = _parse_feature_components(feat)
            if components["G"]:
                g_set.add(components["G"])
            if components["GxG"]:
                gxg_set.add(components["GxG"])
            if components["E"]:
                e_set.add(components["E"])

        results[model] = {
            "G"  : {"count": len(g_set),   "features": g_set},
            "GxG": {"count": len(gxg_set), "features": gxg_set},
            "E"  : {"count": len(e_set),   "features": e_set},
        }

        for ft, data in results[model].items():
            log.info("  %-4s : %d unique features", ft, data["count"])

    return results


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def build_summary_df(results: dict) -> pd.DataFrame:
    """Wide table — one row per model with G_count, GxG_count, E_count, total."""
    rows = []
    for model, ftypes in results.items():
        row   = {"model": model}
        total = 0
        for ft in ("G", "GxG", "E"):
            cnt          = ftypes[ft]["count"]
            row[f"{ft}_count"] = cnt
            total       += cnt
        row["total_count"] = total
        rows.append(row)
    return pd.DataFrame(rows)


def build_long_df(results: dict) -> pd.DataFrame:
    """Long table — one row per (model, feature_type)."""
    rows = []
    for model, ftypes in results.items():
        for ft in ("G", "GxG", "E"):
            rows.append({"model": model, "feature_type": ft, "count": ftypes[ft]["count"]})
    return pd.DataFrame(rows)


def write_feature_lists(results: dict, out_dir: Path) -> None:
    """Write one CSV per (model, feature_type) listing every unique feature."""
    lists_dir = out_dir / "feature_lists"
    lists_dir.mkdir(parents=True, exist_ok=True)

    for model, ftypes in results.items():
        for ft in ("G", "GxG", "E"):
            feats = sorted(ftypes[ft]["features"])
            if feats:
                p = lists_dir / f"{model}_{ft}_features.csv"
                pd.DataFrame(feats, columns=["feature"]).to_csv(p, index=False)
                log.info("  Written: %s  (%d features)", p.name, len(feats))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(pheno_data: str, input_file: str) -> None:
    """Load → filter → parse → count → write outputs."""

    pheno_path  = Path(pheno_data)
    input_path  = (
        Path(input_file) if Path(input_file).is_absolute()
        else pheno_path / input_file
    )
    output_path = pheno_path / "scores"
    output_path.mkdir(parents=True, exist_ok=True)

    log.info("Input  : %s", input_path)
    log.info("Output : %s", output_path)

    # ---- Load --------------------------------------------------------------
    df = pd.read_csv(input_path)
    log.info("Loaded %d rows × %d columns.", *df.shape)

    for col in (MODEL_COL, FEATURE_COL, COEF_COL):
        if col not in df.columns:
            raise ValueError(
                f"Required column '{col}' not found. "
                f"Available: {df.columns.tolist()}"
            )

    df[COEF_COL] = pd.to_numeric(df[COEF_COL], errors="coerce").fillna(0)

    # ---- Count -------------------------------------------------------------
    results = count_features(df)

    # ---- Console summary ---------------------------------------------------
    summary_df = build_summary_df(results)
    print("\n" + "=" * 58)
    print("  PRS Feature Counts  (non-covariate │ non-zero coef)")
    print("=" * 58)
    print(summary_df.to_string(index=False))
    print("=" * 58 + "\n")

    # ---- Write outputs -----------------------------------------------------
    p = output_path / "feature_counts_summary.csv"
    summary_df.to_csv(p, index=False)
    log.info("Summary table  → %s", p)

    p = output_path / "feature_counts_long.csv"
    build_long_df(results).to_csv(p, index=False)
    log.info("Long table     → %s", p)

    write_feature_lists(results, output_path)
    log.info("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Count unique G, GxG, and E features in cardio_product and "
            "cardio_summed PRS models from the combined filtered feature scores CSV."
        )
    )
    parser.add_argument(
        "--pheno_data",
        type=str,
        required=True,
        help="Root phenotype data directory. Outputs go to {pheno_data}/scores/.",
    )
    parser.add_argument(
        "--input_file",
        type=str,
        default="featureScoresReducedFinalModel.combined.filtered.csv",
        help=(
            "Filename or absolute path to the combined filtered CSV. "
            "Resolved relative to --pheno_data if not absolute. "
            "Default: featureScoresReducedFinalModel.combined.filtered.csv"
        ),
    )

    args = parser.parse_args()
    run(pheno_data=args.pheno_data, input_file=args.input_file)
