"""
crosstab_subtype_analysis.py
────────────────────────────
Cross-tabulate genomic (per-cohort) vs clinical (combined) bNMF subtypes.

For each cohort the script:
  1. Loads clinical cluster assignments from
       combinedCohortBnmf_{population}_clinical_only/cluster_assignments.csv
  2. Loads genomic cluster assignments from
       cohortBnmf/{cohort}/[{population}/]cluster_assignments.csv
  3. Merges on IID, optionally filtering to confident assignments.
  4. Produces cross-tab counts, row-%, col-%, chi-square / Cramér's V.
  5. Saves per-cohort CSV tables and a multi-panel heatmap figure.

Usage:
    python scripts/crosstab_subtype_analysis.py \
        --pheno_data results/type2Diabetes/combinedAnalysis \
        --population high_risk \
        --output_dir results/type2Diabetes/combinedAnalysis/scores/crosstab_subtype

    # Or point directly to the scores folder:
    python scripts/crosstab_subtype_analysis.py \
        --scores_dir results/type2Diabetes/combinedAnalysis/scores \
        --population high_risk
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

# ── optional plotting ────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── constants ────────────────────────────────────────────────────────────────
_POPULATION_DIR_NAMES = {"high_risk", "low_control", "both", "separate"}
_SKIP_DIRS = _POPULATION_DIR_NAMES | {"enrichment", "__pycache__"}


# ─────────────────────────────────────────────────────────────────────────────
# Path resolution helpers
# ─────────────────────────────────────────────────────────────────────────────

def _scores_root(pheno_data: str | Path) -> Path:
    """Return the `scores/` directory given a pheno_data path."""
    p = Path(pheno_data)
    for candidate in (p / "scores", p / "combinedAnalysis" / "scores"):
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        f"Cannot find scores/ directory under {pheno_data}. "
        "Pass --scores_dir directly."
    )


def _list_cohort_dirs(cohort_bnmf: Path) -> list[Path]:
    """Return sub-directories of cohortBnmf/ that look like real cohorts."""
    return [
        d for d in sorted(cohort_bnmf.iterdir())
        if d.is_dir() and d.name not in _SKIP_DIRS
    ]


def _resolve_genomic_dir(cohort_dir: Path, population: str) -> Path:
    """
    Return the directory containing feature_loadings.csv / cluster_assignments.csv
    for a given cohort.  Checks `{cohort}/{population}/` first (for separate-pop runs),
    falls back to `{cohort}/`.
    """
    if population not in ("both", "separate"):
        candidate = cohort_dir / population
        if (candidate / "cluster_assignments.csv").exists():
            return candidate
    if (cohort_dir / "cluster_assignments.csv").exists():
        return cohort_dir
    return cohort_dir   # return anyway; caller will raise a better error


def _find_clinical_dir(scores_root: Path, population: str) -> Path:
    """Locate combinedCohortBnmf_{population}_clinical_only/."""
    # Try exact name
    candidate = scores_root / f"combinedCohortBnmf_{population}_clinical_only"
    if candidate.is_dir():
        return candidate
    # Fuzzy: any folder containing 'clinical_only'
    matches = [d for d in scores_root.iterdir()
               if d.is_dir() and "clinical_only" in d.name]
    if len(matches) == 1:
        return matches[0]
    if matches:
        # prefer the one matching population
        for m in matches:
            if population in m.name:
                return m
        return matches[0]
    raise FileNotFoundError(
        f"No clinical-only NMF directory found in {scores_root}. "
        f"Expected: combinedCohortBnmf_{population}_clinical_only"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Load helpers
# ─────────────────────────────────────────────────────────────────────────────

def _load_assignments(path: Path, id_col: str | None = "IID") -> pd.DataFrame:
    """
    Load cluster_assignments.csv.
    Handles both:
      - named first column (id_col='IID')
      - unnamed index (id_col=None → use index)
    Returns DataFrame with 'IID', 'hard_cluster', 'cluster_label',
    'membership_entropy', and all soft-membership columns.
    """
    df = pd.read_csv(path)
    if id_col and id_col in df.columns:
        df = df.rename(columns={id_col: "IID"})
    elif df.columns[0] in ("", "Unnamed: 0"):
        df = df.rename(columns={df.columns[0]: "IID"})
    else:
        # first column might already be unnamed index
        df.insert(0, "IID", df.index)
        df = df.reset_index(drop=True)
    df["IID"] = df["IID"].astype(str)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Statistics
# ─────────────────────────────────────────────────────────────────────────────

def cramers_v(contingency: pd.DataFrame) -> float:
    """Cramér's V effect-size for a contingency table."""
    chi2, _, _, _ = chi2_contingency(contingency, correction=False)
    n = contingency.values.sum()
    r, c = contingency.shape
    return float(np.sqrt(chi2 / (n * (min(r, c) - 1)))) if min(r, c) > 1 else 0.0


def run_chi2(contingency: pd.DataFrame) -> dict:
    """Return chi2, p, dof, cramers_v for a contingency table."""
    chi2, p, dof, _ = chi2_contingency(contingency, correction=False)
    return {
        "chi2": round(chi2, 4),
        "p_value": float(f"{p:.4e}"),
        "dof": int(dof),
        "cramers_v": round(cramers_v(contingency), 4),
        "n": int(contingency.values.sum()),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Soft-membership overlap
# ─────────────────────────────────────────────────────────────────────────────

def _membership_cols(df: pd.DataFrame, prefix: str = "membership_") -> list[str]:
    return [c for c in df.columns if c.startswith(prefix)]


def soft_crosstab(
    gen_df: pd.DataFrame,
    clin_df: pd.DataFrame,
    gen_label_col: str = "cluster_label",
    clin_label_col: str = "cluster_label",
) -> pd.DataFrame:
    """
    Compute a continuous co-membership matrix: for each pair (genomic cluster i,
    clinical cluster j), sum over all shared samples the product
    w_gen[i] × w_clin[j].  Gives fractional counts proportional to how strongly
    each sample belongs to each cluster pair.
    """
    merged = pd.merge(
        gen_df, clin_df, on="IID", suffixes=("_gen", "_clin")
    )
    gen_mem = _membership_cols(merged, "membership_")
    # Disambiguate when both have membership_ cols after merge
    gen_mem_cols = [c for c in merged.columns if c.startswith("membership_") and c.endswith("_gen")]
    clin_mem_cols = [c for c in merged.columns if c.startswith("membership_") and c.endswith("_clin")]

    if not gen_mem_cols or not clin_mem_cols:
        return pd.DataFrame()

    gen_labels = [c.replace("membership_", "").replace("_gen", "") for c in gen_mem_cols]
    clin_labels = [c.replace("membership_", "").replace("_clin", "") for c in clin_mem_cols]

    mat = np.zeros((len(gen_labels), len(clin_labels)))
    for gi, gc in enumerate(gen_mem_cols):
        for ci, cc in enumerate(clin_mem_cols):
            mat[gi, ci] = (merged[gc] * merged[cc]).sum()

    return pd.DataFrame(mat, index=gen_labels, columns=clin_labels)


# ─────────────────────────────────────────────────────────────────────────────
# Per-cohort cross-tabulation
# ─────────────────────────────────────────────────────────────────────────────

def crosstab_cohort(
    cohort: str,
    genomic_dir: Path,
    clin_df: pd.DataFrame,
    population: str,
    min_entropy: float | None = None,
    max_entropy: float | None = None,
) -> dict | None:
    """
    Build cross-tab for one cohort.  Returns a dict with keys:
      counts, rowpct, colpct, soft, stats, n_merged, n_filtered, cohort
    Returns None if cluster_assignments.csv is missing.
    """
    asgn_path = genomic_dir / "cluster_assignments.csv"
    if not asgn_path.exists():
        print(f"  [SKIP] {cohort}: no cluster_assignments.csv at {asgn_path}", flush=True)
        return None

    gen_df = _load_assignments(asgn_path, id_col="IID")
    n_raw = len(gen_df)

    # Optional entropy filter (exclude ambiguous assignments)
    if "membership_entropy" in gen_df.columns:
        if min_entropy is not None:
            gen_df = gen_df[gen_df["membership_entropy"] >= min_entropy]
        if max_entropy is not None:
            gen_df = gen_df[gen_df["membership_entropy"] <= max_entropy]
    n_filtered = n_raw - len(gen_df)

    # Merge with clinical
    merged = pd.merge(
        gen_df[["IID", "hard_cluster", "cluster_label", "membership_entropy"]
               + _membership_cols(gen_df)].copy(),
        clin_df[["IID", "hard_cluster", "cluster_label"]
                + _membership_cols(clin_df)].rename(
            columns={"hard_cluster": "clin_hard_cluster",
                     "cluster_label": "clin_cluster_label"}
        ),
        on="IID",
        how="inner",
    )

    if merged.empty:
        print(f"  [SKIP] {cohort}: no IID overlap after merge", flush=True)
        return None

    # Hard-assignment cross-tab
    gen_col = "cluster_label"
    clin_col = "clin_cluster_label"
    counts = pd.crosstab(
        merged[gen_col].rename("genomic_cluster"),
        merged[clin_col].rename("clinical_cluster"),
    )
    counts.index.name = "genomic_cluster"
    counts.columns.name = "clinical_cluster"

    rowpct = counts.div(counts.sum(axis=1), axis=0).mul(100).round(1)
    colpct = counts.div(counts.sum(axis=0), axis=1).mul(100).round(1)

    # Soft membership cross-tab
    gen_mem_df = gen_df.copy()
    gen_mem_df = gen_mem_df.rename(
        columns={c: c for c in gen_mem_df.columns}  # keep as-is, suffix added in soft_crosstab
    )
    clin_mem_df = clin_df.copy()
    soft = soft_crosstab(gen_mem_df, clin_mem_df)

    # Statistics
    stats = run_chi2(counts)
    stats["cohort"] = cohort
    stats["n_merged"] = len(merged)
    stats["n_filtered"] = n_filtered
    stats["genomic_k"] = counts.shape[0]
    stats["clinical_k"] = counts.shape[1]

    return {
        "cohort": cohort,
        "counts": counts,
        "rowpct": rowpct,
        "colpct": colpct,
        "soft": soft,
        "stats": stats,
        "merged": merged,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def _heatmap_panel(
    ax,
    data: pd.DataFrame,
    title: str,
    fmt: str = ".0f",
    cmap: str = "Blues",
    vmin=None,
    vmax=None,
    annot_size: int = 8,
):
    """Draw a single annotated heatmap on ax."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    arr = data.values.astype(float)
    im = ax.imshow(arr, cmap=cmap, aspect="auto",
                   norm=Normalize(vmin=vmin or arr.min(), vmax=vmax or arr.max()))
    ax.set_xticks(range(data.shape[1]))
    ax.set_yticks(range(data.shape[0]))
    ax.set_xticklabels(data.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(data.index, fontsize=7)
    ax.set_xlabel("Clinical cluster", fontsize=8)
    ax.set_ylabel("Genomic cluster", fontsize=8)
    ax.set_title(title, fontsize=9, pad=4)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            txt = format(val, fmt)
            color = "white" if val > (arr.max() * 0.6) else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=annot_size, color=color)
    return im


def make_heatmap_figure(results: list[dict], output_path: Path):
    """
    Multi-panel figure: one row per cohort, three columns:
      counts | row% | col%
    """
    if not HAS_MPL:
        print("[WARN] matplotlib not available — skipping heatmap figure")
        return

    n_cohorts = len(results)
    fig, axes = plt.subplots(
        n_cohorts, 3,
        figsize=(14, 4 * n_cohorts),
        squeeze=False,
    )
    for row_idx, res in enumerate(results):
        cohort = res["cohort"]
        stats = res["stats"]
        subtitle = (
            f"{cohort}  (n={stats['n_merged']}, "
            f"χ²={stats['chi2']}, p={stats['p_value']:.2e}, "
            f"V={stats['cramers_v']})"
        )
        _heatmap_panel(axes[row_idx, 0], res["counts"],
                       f"Counts — {subtitle}", fmt=".0f", cmap="Blues")
        _heatmap_panel(axes[row_idx, 1], res["rowpct"],
                       "Row % (genomic→clinical)", fmt=".1f", cmap="YlOrRd",
                       vmin=0, vmax=100)
        _heatmap_panel(axes[row_idx, 2], res["colpct"],
                       "Col % (clinical→genomic)", fmt=".1f", cmap="YlGn",
                       vmin=0, vmax=100)

    fig.suptitle("Genomic × Clinical Subtype Cross-Tabulation", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved heatmap → {output_path}")


def make_soft_figure(results: list[dict], output_path: Path):
    """Soft-membership co-occurrence heatmap (one panel per cohort)."""
    if not HAS_MPL:
        return
    valid = [r for r in results if r["soft"] is not None and not r["soft"].empty]
    if not valid:
        return

    n = len(valid)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), squeeze=False)
    for idx, res in enumerate(valid):
        soft_norm = res["soft"].div(res["soft"].sum().sum())
        _heatmap_panel(axes[0, idx], soft_norm * 100,
                       f"Soft co-membership %\n{res['cohort']}", fmt=".1f",
                       cmap="Purples", vmin=0, vmax=None)
    fig.suptitle("Soft Membership Co-occurrence (genomic × clinical)", fontsize=11)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved soft heatmap → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Cross-tabulate genomic (per-cohort) vs clinical (combined) bNMF subtypes."
    )
    loc = parser.add_argument_group("Paths")
    loc.add_argument(
        "--pheno_data", default=None,
        help="Path to phenotype data root (e.g. results/type2Diabetes/combinedAnalysis). "
             "Used to auto-detect scores_dir.",
    )
    loc.add_argument(
        "--scores_dir", default=None,
        help="Path to scores/ directory (overrides --pheno_data auto-detection).",
    )
    loc.add_argument(
        "--output_dir", default=None,
        help="Where to save outputs. Defaults to {scores_dir}/crosstab_subtype_{population}.",
    )
    loc.add_argument(
        "--cohorts", nargs="+", default=None,
        help="Specific cohort subdirectories to process. Default: all found in cohortBnmf/.",
    )

    mode = parser.add_argument_group("Analysis mode")
    mode.add_argument(
        "--population", default="high_risk",
        choices=["high_risk", "low_control", "both"],
        help="Population used for both genomic and clinical runs (default: high_risk).",
    )
    mode.add_argument(
        "--min_entropy", type=float, default=None,
        help="Exclude samples with membership_entropy < this value (low confidence).",
    )
    mode.add_argument(
        "--max_entropy", type=float, default=None,
        help="Exclude samples with membership_entropy > this value (high ambiguity).",
    )

    parser.add_argument(
        "--no_plots", action="store_true",
        help="Skip heatmap figures.",
    )

    args = parser.parse_args(argv)

    # ── Resolve scores_dir ────────────────────────────────────────────────
    if args.scores_dir:
        scores_root = Path(args.scores_dir)
    elif args.pheno_data:
        scores_root = _scores_root(args.pheno_data)
    else:
        parser.error("Provide --pheno_data or --scores_dir.")

    print(f"Scores root: {scores_root}")

    cohort_bnmf = scores_root / "cohortBnmf"
    if not cohort_bnmf.is_dir():
        sys.exit(f"[ERROR] cohortBnmf/ not found at {cohort_bnmf}")

    # ── Clinical assignments ───────────────────────────────────────────────
    clin_dir = _find_clinical_dir(scores_root, args.population)
    clin_asgn = clin_dir / "cluster_assignments.csv"
    if not clin_asgn.exists():
        sys.exit(f"[ERROR] Clinical assignments not found: {clin_asgn}")

    print(f"Clinical NMF dir: {clin_dir.name}")
    clin_df = _load_assignments(clin_asgn, id_col=None)
    print(f"  {len(clin_df)} clinical assignments loaded")

    # ── Output dir ────────────────────────────────────────────────────────
    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        out_dir = scores_root / f"crosstab_subtype_{args.population}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output dir: {out_dir}")

    # ── Cohorts to process ───────────────────────────────────────────────
    cohort_dirs = _list_cohort_dirs(cohort_bnmf)
    if args.cohorts:
        cohort_dirs = [d for d in cohort_dirs if d.name in args.cohorts]
    print(f"\nProcessing {len(cohort_dirs)} cohort(s): "
          f"{[d.name for d in cohort_dirs]}")

    # ── Run per-cohort ────────────────────────────────────────────────────
    all_results = []
    summary_rows = []

    for cohort_dir in cohort_dirs:
        cohort = cohort_dir.name
        print(f"\n── {cohort} ──────────────────────────────────", flush=True)
        genomic_dir = _resolve_genomic_dir(cohort_dir, args.population)
        res = crosstab_cohort(
            cohort=cohort,
            genomic_dir=genomic_dir,
            clin_df=clin_df,
            population=args.population,
            min_entropy=args.min_entropy,
            max_entropy=args.max_entropy,
        )
        if res is None:
            continue

        all_results.append(res)
        stats = res["stats"]
        summary_rows.append(stats)

        # Print counts table
        print(f"  n={stats['n_merged']}  "
              f"χ²={stats['chi2']}  p={stats['p_value']:.3e}  "
              f"V={stats['cramers_v']}")
        print(res["counts"].to_string())

        # Save per-cohort tables
        prefix = out_dir / f"crosstab_{cohort}"
        res["counts"].to_csv(f"{prefix}_counts.csv")
        res["rowpct"].to_csv(f"{prefix}_rowpct.csv")
        res["colpct"].to_csv(f"{prefix}_colpct.csv")
        if not res["soft"].empty:
            res["soft"].to_csv(f"{prefix}_soft.csv")

        # Save merged IID-level data
        # membership_entropy may be suffixed _x after merge
        entropy_col = next(
            (c for c in ("membership_entropy", "membership_entropy_x") if c in res["merged"].columns),
            None,
        )
        save_cols = ["IID", "cluster_label", "clin_cluster_label"]
        if entropy_col:
            save_cols.append(entropy_col)
        res["merged"][save_cols].rename(
            columns={entropy_col: "membership_entropy"} if entropy_col else {}
        ).to_csv(f"{prefix}_iid_assignments.csv", index=False)

    if not all_results:
        print("\n[WARN] No results produced — check cohort paths and IID overlap.")
        return

    # ── Summary table ─────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows).set_index("cohort")
    col_order = ["n_merged", "n_filtered", "genomic_k", "clinical_k",
                 "chi2", "p_value", "dof", "cramers_v"]
    summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]
    summary_df.to_csv(out_dir / "crosstab_summary.csv")
    print(f"\n{'─'*60}")
    print("Summary:")
    print(summary_df.to_string())
    print(f"\nSaved summary → {out_dir / 'crosstab_summary.csv'}")

    # ── Figures ───────────────────────────────────────────────────────────
    if not args.no_plots:
        make_heatmap_figure(all_results, out_dir / "crosstab_heatmap.pdf")
        make_soft_figure(all_results, out_dir / "crosstab_soft_heatmap.pdf")

    print(f"\nDone. Output in: {out_dir}")


if __name__ == "__main__":
    main()
