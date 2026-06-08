#!/usr/bin/env python3
"""
cluster_subtype_dotmatrix.py
============================
Bubble (dot-matrix) chart: functional clusters × T2D subtypes,
with bubble size ∝ % of cohort features and colour per PRS model.

Reads:
  results/type2Diabetes/combinedAnalysis/scores/functionalEnrichment/
      cluster_subtype_cohort_weighted.csv

Usage
-----
python scripts/cluster_subtype_dotmatrix.py \\
    [--data   PATH_TO_CSV] \\
    [--out_dir OUTPUT_DIR]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd

# ── Cohort → PRS model key ────────────────────────────────────────────────────
COHORT_MODEL = {
    'main':           'PRS_G',
    'epi_summed':     'PRS_G+G',
    'epi_product':    'PRS_GxG',
    'cardio_summed':  'PRS_G+GxE',
    'cardio_product': 'PRS_GxGxE',
}

# ── Cluster key → display name (merges beta_cell_1+2 and body_fat+obesity) ───
CLUSTER_MERGE = {
    'hla_autoimmune':     'HLA autoimmune',
    'beta_cell_1':        'Beta-cell (1+2)',
    'beta_cell_2':        'Beta-cell (1+2)',
    'proinsulin':         'Proinsulin',
    'res_glycaemic':      'Residual glycaemia',
    'hyper_insulin':      'Hyper-ins. secretion',
    'liver_lipid':        'Liver / lipid metab.',
    'body_fat':           'Obesity / body fat',
    'obesity':            'Obesity / body fat',
    'metabolic_syndrome': 'Metabolic syndrome',
    'lipodystrophy':      'Lipodystrophy',
}

# ── Biological-pathway ordering (top → bottom) ────────────────────────────────
CLUSTERS = [
    'HLA autoimmune',
    'Beta-cell (1+2)',
    'Proinsulin',
    'Residual glycaemia',
    'Hyper-ins. secretion',
    'Liver / lipid metab.',
    'Obesity / body fat',
    'Metabolic syndrome',
    'Lipodystrophy',
]

SUBTYPES = ['SAID', 'SIDD', 'MOD', 'MARD', 'SIRD']

MODELS = ['PRS_G', 'PRS_G+G', 'PRS_GxG', 'PRS_G+GxE', 'PRS_GxGxE']
MODEL_LABELS = {
    'PRS_G':     'PRS₂G  (additive)',
    'PRS_G+G':   'PRS₂G+G (G burden)',
    'PRS_GxG':   'PRS₂GxG (epistasis)',
    'PRS_G+GxE': 'PRS₂G+GxE',
    'PRS_GxGxE': 'PRS₂GxGxE',
}
MODEL_COLORS = {
    'PRS_G':     '#E69F00',  # amber        (main)
    'PRS_G+G':   '#0073A8',  # dark blue    (epi_summed)
    'PRS_GxG':   '#56B4E9',  # sky blue     (epi_product)
    'PRS_G+GxE': '#009E73',  # teal-green   (cardio_summed)
    'PRS_GxGxE': '#006B52',  # forest green (cardio_product)
}

# Background colour per section
CLUSTER_BG = {
    'HLA autoimmune':       '#fde8e8',  # soft pink   (immune)
    'Beta-cell (1+2)':      '#e4eff9',  # soft blue
    'Proinsulin':           '#e4eff9',
    'Residual glycaemia':   '#e4eff9',
    'Hyper-ins. secretion': '#e4eff9',
    'Liver / lipid metab.': '#e8f5e3',  # soft green  (metabolic)
    'Obesity / body fat':   '#e8f5e3',
    'Metabolic syndrome':   '#e8f5e3',
    'Lipodystrophy':        '#e8f5e3',
}

# Area (pt²) for 50 % weight — sets the overall bubble scale
MAX_AREA = 900


def pct_to_area(pct: float) -> float:
    return pct / 50.0 * MAX_AREA


def load_data(csv_path: str) -> pd.DataFrame:
    raw = pd.read_csv(csv_path)
    raw = raw[raw['cohort'].isin(COHORT_MODEL)]
    raw = raw[raw['cluster_key'].isin(CLUSTER_MERGE)].copy()

    raw['model']   = raw['cohort'].map(COHORT_MODEL)
    raw['cluster'] = raw['cluster_key'].map(CLUSTER_MERGE)

    grouped = (
        raw.groupby(['model', 'cluster', 'subtype'])
           .agg(weighted_n=('weighted_n', 'sum'),
                cohort_total_w=('cohort_total_w', 'first'))
           .reset_index()
    )
    grouped['weight'] = (
        grouped['weighted_n'] / grouped['cohort_total_w'] * 100
    ).round(1)
    return grouped[grouped['weight'] > 0].copy()


def make_figure(df: pd.DataFrame, out_dir: str) -> None:
    n_clusters = len(CLUSTERS)
    n_subtypes = len(SUBTYPES)

    cluster_idx = {c: i for i, c in enumerate(CLUSTERS)}
    subtype_idx = {s: i for i, s in enumerate(SUBTYPES)}
    model_idx   = {m: i for i, m in enumerate(MODELS)}

    # Horizontal jitter within each subtype column (5 models side-by-side)
    offsets = np.linspace(-0.22, 0.22, len(MODELS))

    fig, ax = plt.subplots(figsize=(13.5, 7.5))

    # ── background row shading ────────────────────────────────────────────────
    for ci, cname in enumerate(CLUSTERS):
        y = n_clusters - 1 - ci
        ax.axhspan(y - 0.5, y + 0.5,
                   color=CLUSTER_BG[cname], zorder=0, linewidth=0)

    # ── subtle horizontal grid lines ──────────────────────────────────────────
    for ci in range(n_clusters):
        ax.axhline(n_clusters - 1 - ci, color='white', lw=1.2, zorder=1)

    # ── bubbles ───────────────────────────────────────────────────────────────
    for _, row in df.iterrows():
        if row['cluster'] not in cluster_idx or row['subtype'] not in subtype_idx:
            continue
        x = subtype_idx[row['subtype']] + offsets[model_idx[row['model']]]
        y = n_clusters - 1 - cluster_idx[row['cluster']]
        ax.scatter(x, y,
                   s=pct_to_area(row['weight']),
                   color=MODEL_COLORS[row['model']],
                   alpha=0.85,
                   edgecolors='none',
                   zorder=4)

    # ── axes ──────────────────────────────────────────────────────────────────
    ax.set_xticks(range(n_subtypes))
    ax.set_xticklabels(SUBTYPES, fontsize=13, fontweight='bold')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(list(reversed(CLUSTERS)), fontsize=10.5)

    ax.set_xlim(-0.6, n_subtypes - 0.4)
    ax.set_ylim(-0.55, n_clusters - 0.45)
    ax.tick_params(bottom=False, top=False, left=False)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # ── section dashed separators ─────────────────────────────────────────────
    # Between HLA autoimmune (ci=0, y=8) and beta-cell block
    immune_y    = n_clusters - 1 - 0   # y = 8
    # Between Hyper-ins. secretion (ci=4, y=4) and metabolic block
    betacell_y  = n_clusters - 1 - 4   # y = 4

    for y_sep in (immune_y - 0.5, betacell_y - 0.5):
        ax.axhline(y_sep, color='#aaaaaa', lw=1.0, ls='--', zorder=2)

    # ── section italic labels (left margin, inside axes) ──────────────────────
    ax.text(-0.64, (immune_y - 1 + betacell_y) / 2,
            'β-cell / glycaemia',
            va='center', ha='right', fontsize=9, color='#2477b5', style='italic',
            clip_on=False)
    ax.text(-0.64, (betacell_y - 1) / 2,
            'Metabolic',
            va='center', ha='right', fontsize=9, color='#2a7a2a', style='italic',
            clip_on=False)

    # ── title ─────────────────────────────────────────────────────────────────
    ax.set_title(
        'Functional cluster → T2D subtype: weighted contributions by PRS model',
        fontsize=12, pad=18,
    )

    # ── model legend (upper right, outside axes) ──────────────────────────────
    model_handles = [
        mpatches.Patch(color=MODEL_COLORS[m], label=MODEL_LABELS[m])
        for m in MODELS
    ]
    leg_model = fig.legend(
        handles=model_handles,
        title='Model',
        title_fontsize=9.5,
        loc='upper right',
        bbox_to_anchor=(1.00, 0.97),
        frameon=True,
        fontsize=9,
        edgecolor='#cccccc',
        handlelength=1.2,
    )

    # ── weight size legend (lower right, outside axes) ────────────────────────
    size_weights = [5, 15, 30, 50]
    size_handles = [
        mlines.Line2D([], [], marker='o', color='#888888', alpha=0.80,
                      markersize=np.sqrt(pct_to_area(w) / np.pi),
                      linestyle='None', label=f'{w}%')
        for w in size_weights
    ]
    fig.legend(
        handles=size_handles,
        title='Weight (%)',
        title_fontsize=9.5,
        loc='lower right',
        bbox_to_anchor=(1.00, 0.08),
        frameon=True,
        fontsize=9,
        edgecolor='#cccccc',
        handletextpad=1.5,
    )

    # ── caption ───────────────────────────────────────────────────────────────
    fig.text(
        0.01, -0.03,
        ('Bubble area ∝ weighted % of model features attributed to that cluster–subtype pair.\n'
         'PRS₂G+GxE distributes evenly across all subtypes; '
         'PRS₂GxGxE is maximally polarised (HLA→SAID, lipodystrophy→SIRD).'),
        fontsize=7.5, color='#555555', va='top',
    )

    fig.subplots_adjust(left=0.18, right=0.68, top=0.88, bottom=0.06)

    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        out = os.path.join(out_dir, f'cluster_subtype_dotmatrix.{ext}')
        fig.savefig(out, dpi=180, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--data',
        default=(
            'results/type2Diabetes/combinedAnalysis/scores'
            '/functionalEnrichment/cluster_subtype_cohort_weighted.csv'
        ),
        help='Path to cluster_subtype_cohort_weighted.csv',
    )
    parser.add_argument(
        '--out_dir',
        default='results/type2Diabetes/subtypeResearch/enrichment/plots',
        help='Output directory for PNG and PDF',
    )
    args = parser.parse_args()
    df = load_data(args.data)
    make_figure(df, args.out_dir)


if __name__ == '__main__':
    main()
