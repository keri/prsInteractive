"""
PRS Model Comparison: Box Plots + Odds Ratio / Prevalence across Deciles
Input : combinedPRSGroups.holdout.filtered.csv
Output: {pheno_data}/figures/modelComparisons/
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator
from scipy.stats import mannwhitneyu
from matplotlib.colors import to_rgba
import warnings
warnings.filterwarnings("ignore")

# ─── CONFIG ──────────────────────────────────────────────────────────────────
# Set PHENO_DATA env-var or edit the default here
PHENO_DATA   = os.environ.get("PHENO_DATA", "/path/to/pheno_data")
INPUT_FILE   = os.path.join(PHENO_DATA, "scores","combinedPRSGroups.holdout.filtered.csv")
OUTPUT_DIR   = os.path.join(PHENO_DATA, "figures", "modelComparisons","CompositeComparison")
N_BINS       = 1000        # number of bins in bin_combined / bin_main
N_DECILES    = 10          # deciles for OR / prevalence
FIGSIZE_BOX  = (14, 6)
FIGSIZE_OR   = (14, 6)
DPI          = 150

# ─── COLOUR PALETTE ──────────────────────────────────────────────────────────
COHORT_COLORS_EXTENDED = {
  'main': '#E69F00',      # Orange
  'epi': '#56B4E9',       # Sky blue
  'epi+main': '#CC79A7',  # Pinkish purple
  'cardio': '#009E73',    # Bluish green
  'all': '#F0E442',       # Yellow
  'combined': '#D55E00',  # Vermillion
  # Product variants (same as base)
  'epi_product': '#56B4E9',
  'cardio_product': '#009E73',
  'epi+main_product': '#CC79A7',
  'all_product': '#F0E442',
  'main_product': '#E69F00',
  # Summed variants (darker shades)
  'epi_summed': '#0073A8',      # Darker blue
  'cardio_summed': '#006B52',   # Darker green
  'epi+main_summed': '#9F5580', # Darker purple
  'all_summed': '#C7B800',      # Darker yellow
  'main_summed': '#C77D00',     # Darker orange
}


COLOR_COMBINED = COHORT_COLORS_EXTENDED["combined"]
COLOR_MAIN     = COHORT_COLORS_EXTENDED["main"]

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def make_output_dir(path: str):
    os.makedirs(path, exist_ok=True)


def detect_status_column(df: pd.DataFrame) -> str:
    """Try common names for the binary case/control column."""
    candidates = ["status", "phenotype", "case", "label",
                  "affected", "disease", "y", "outcome"]
    cols_lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in cols_lower:
            return cols_lower[c]
    # fall back: first binary column that isn't a bin
    for col in df.columns:
        if col.lower().startswith("bin"):
            continue
        vals = df[col].dropna().unique()
        if set(vals).issubset({0, 1, 0.0, 1.0}):
            return col
    raise ValueError(
        "Cannot find a binary case/control column. "
        "Set it manually with STATUS_COL variable."
    )


def bin_to_score(series: pd.Series, n_bins: int = N_BINS) -> pd.Series:
    """Normalise 1-based bin number to 0-1 continuous PRS-like score."""
    return (series - 1) / (n_bins - 1)


def decile_groups(series: pd.Series, n: int = N_DECILES) -> pd.Series:
    """Assign each row to a 1-indexed decile based on a numeric series."""
    return pd.qcut(series, q=n, labels=range(1, n + 1)).astype(int)


def odds_ratio_per_group(df: pd.DataFrame,
                         group_col: str,
                         status_col: str,
                         reference: int = 1):
    """
    Compute OR and 95 % CI (Woolf / log method) for each group vs. reference.
    Returns a DataFrame with columns: group, OR, lower, upper, prevalence.
    """
    records = []
    ref = df[df[group_col] == reference]
    a0 = (ref[status_col] == 2).sum()
    b0 = (ref[status_col] == 1).sum()

    for g in sorted(df[group_col].unique()):
        grp   = df[df[group_col] == g]
        cases = (grp[status_col] == 2).sum()
        ctrls = (grp[status_col] == 1).sum()
        prev  = cases / len(grp) if len(grp) > 0 else np.nan

        if g == reference:
            records.append(dict(group=g, OR=1.0, lower=1.0,
                                upper=1.0, prevalence=prev))
            continue

        # Haldane–Anscombe correction when zeros present
        a1, b1, a0c, b0c = (cases + 0.5, ctrls + 0.5,
                             a0 + 0.5, b0 + 0.5)
        if a0 == 0 or b0 == 0 or cases == 0 or ctrls == 0:
          a1, b1, a0c, b0c = cases + 0.5, ctrls + 0.5, a0 + 0.5, b0 + 0.5
        else:
          a1, b1, a0c, b0c = cases, ctrls, a0, b0

        log_or = np.log(a1 * b0c) - np.log(b1 * a0c)
        se_log = np.sqrt(1/a1 + 1/b1 + 1/a0c + 1/b0c)
        OR     = np.exp(log_or)
        lo     = np.exp(log_or - 1.96 * se_log)
        hi     = np.exp(log_or + 1.96 * se_log)
        records.append(dict(group=g, OR=OR, lower=lo,
                            upper=hi, prevalence=prev))

    return pd.DataFrame(records)


# ─── PLOTTING ────────────────────────────────────────────────────────────────
def plot_boxplots(df: pd.DataFrame, status_col: str, out_dir: str):
    """
    Side-by-side box plots for bin_combined and bin_main,
    split by case (1) / control (0).
    """
    score_cols = {
        "bin_combined": ("combined", COLOR_COMBINED),
        "bin_main"    : ("main",     COLOR_MAIN),
    }
    # convert bins → 0–1 score for comparability
    for col in score_cols:
        df[col + "_score"] = bin_to_score(df[col])

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE_BOX, sharey=True)
    fig.suptitle("PRS Score Distribution: Cases vs Controls",
                 fontsize=14, fontweight="bold", y=1.01)

    for ax, (col, (label, color)) in zip(axes, score_cols.items()):
        score_col = col + "_score"
        cases_data = df.loc[df[status_col] == 2, score_col].dropna()
        ctrl_data  = df.loc[df[status_col] == 1, score_col].dropna()

        data    = [ctrl_data, cases_data]
        labels  = ["Controls\n(0)", "Cases\n(1)"]
        colors = [
          (*to_rgba(color)[:3], 0.30),   # controls — 30% opacity
          (*to_rgba(color)[:3], 0.90),   # cases    — 90% opacity
        ]

        bp = ax.boxplot(
            data,
            patch_artist=True,
            notch=False,
            widths=0.45,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            flierprops=dict(marker="o", markersize=2,
                            alpha=0.3, markeredgewidth=0),
        )
        for patch, clr in zip(bp["boxes"], colors):
            patch.set_facecolor(clr)


        # Mann-Whitney U
        stat, p = mannwhitneyu(cases_data, ctrl_data, alternative="two-sided")
        p_str = f"p = {p:.2e}" if p >= 1e-300 else "p < 1e-300"
        ax.set_title(f"{col}\n{p_str}", fontsize=11)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel("Normalised PRS Score (0 – 1)", fontsize=10)
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        ax.spines[["top", "right"]].set_visible(False)

        # legend patches
        patches = [
            mpatches.Patch(facecolor=color + "55", label="Controls"),
            mpatches.Patch(facecolor=color,         label="Cases"),
        ]
        ax.legend(handles=patches, fontsize=9, loc="upper left")

    plt.tight_layout()
    out_path = os.path.join(out_dir, "boxplot_prs_cases_controls.png")
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {out_path}")


def plot_or_and_prevalence(df: pd.DataFrame, status_col: str, out_dir: str):
    """
    Odds ratio and prevalence over 10 deciles for bin_combined and bin_main.
    Two-panel figure: top = OR, bottom = prevalence.
    """
    models = [
        ("bin_combined", "Combined Model", COLOR_COMBINED),
        ("bin_main",     "Main Model",     COLOR_MAIN),
    ]

    fig, (ax_or, ax_prev) = plt.subplots(
        2, 1, figsize=FIGSIZE_OR, sharex=True,
        gridspec_kw={"hspace": 0.05}
    )
    fig.suptitle("Odds Ratio & Prevalence across PRS Deciles",
                 fontsize=14, fontweight="bold")

    x = np.arange(1, N_DECILES + 1)

    for col, label, color in models:
        df[f"{col}_decile"] = decile_groups(df[col])
        stats = odds_ratio_per_group(df, f"{col}_decile", status_col)

        # ── OR panel ─────────────────────────────────────────────────────
        ax_or.plot(stats["group"], stats["OR"],
                   marker="o", color=color, label=label, linewidth=2)
        ax_or.fill_between(
            stats["group"], stats["lower"], stats["upper"],
            color=color, alpha=0.15
        )

        # ── Prevalence panel ─────────────────────────────────────────────
        ax_prev.plot(stats["group"], stats["prevalence"] * 100,
                     marker="s", color=color, label=label,
                     linewidth=2, linestyle="--")

    # OR reference line
    ax_or.axhline(1, color="grey", linewidth=1, linestyle=":")
    ax_or.set_ylabel("Odds Ratio vs Decile 1", fontsize=10)
    ax_or.set_yscale("log")
    ax_or.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda y, _: f"{y:g}")
    )
    ax_or.legend(fontsize=9, loc="upper left")
    ax_or.grid(axis="y", alpha=0.25, linestyle="--")
    ax_or.spines[["top", "right"]].set_visible(False)

    ax_prev.set_xlabel("PRS Decile", fontsize=10)
    ax_prev.set_ylabel("Prevalence (%)", fontsize=10)
    ax_prev.set_xticks(x)
    ax_prev.legend(fontsize=9, loc="upper left")
    ax_prev.grid(axis="y", alpha=0.25, linestyle="--")
    ax_prev.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "or_prevalence_deciles.png")
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {out_path}")


def save_decile_stats(df: pd.DataFrame, status_col: str, out_dir: str):
    """Export decile-level OR / prevalence tables to CSV."""
    for col, label in [("bin_combined", "combined"), ("bin_main", "main")]:
        col_d  = f"{col}_decile"
        if col_d not in df.columns:
            df[col_d] = decile_groups(df[col])
        stats = odds_ratio_per_group(df, col_d, status_col)
        stats.columns = ["decile", "OR", "CI_lower", "CI_upper", "prevalence"]
        out_path = os.path.join(out_dir, f"decile_stats_{label}.csv")
        stats.to_csv(out_path, index=False, float_format="%.6f")
        print(f"  [saved] {out_path}")


# ─── MAIN ────────────────────────────────────────────────────────────────────
def main():
    make_output_dir(OUTPUT_DIR)

    print(f"Loading: {INPUT_FILE}")
    if not os.path.isfile(INPUT_FILE):
        sys.exit(
            f"ERROR: File not found: {INPUT_FILE}\n"
            "Set the PHENO_DATA environment variable to the correct directory."
        )

    df = pd.read_csv(INPUT_FILE)
    print(f"  Rows: {len(df):,}   Columns: {list(df.columns)}")

    # ── Validate required columns ─────────────────────────────────────────
    for col in ("bin_combined", "bin_main"):
        if col not in df.columns:
            sys.exit(f"ERROR: Expected column '{col}' not found in CSV.")

    status_col = detect_status_column(df)
    print(f"  Status column: '{status_col}'  "
          f"(cases={( df[status_col]==1).sum():,}  "
          f"controls={(df[status_col]==0).sum():,})")

    # Drop rows missing essential columns
    df = df.dropna(subset=["bin_combined", "bin_main", status_col])
    df[status_col] = df[status_col].astype(int)

    # ── Plots ─────────────────────────────────────────────────────────────
    print("\nGenerating box plots …")
    plot_boxplots(df, status_col, OUTPUT_DIR)

    print("Generating OR / prevalence plots …")
    plot_or_and_prevalence(df, status_col, OUTPUT_DIR)

    print("Saving decile statistics …")
    save_decile_stats(df, status_col, OUTPUT_DIR)

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
  