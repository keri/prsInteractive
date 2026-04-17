#!/usr/bin/env python3
"""
cluster_subtype_dotmatrix.py
============================
Bubble (dot-matrix) chart: functional clusters × T2D subtypes,
with bubble size ∝ % of cohort features and colour per PRS model.

Reads from:
  results/type2Diabetes/subtypeResearch/enrichment/cluster_subtype_cohort_weighted.csv

Usage
-----
python scripts/cluster_subtype_dotmatrix.py \
    [--data   PATH_TO_CSV] \
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

# ── Cohort → PRS model label ──────────────────────────────────────────────────
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
    'PRS_G':       'PRS\u2082G  (additive)',
    'PRS_G+G':     'PRS\u2082G+G  (G burden)',
    'PRS_GxG':     'PRS\u2082GxG  (epistasis)',
    'PRS_G+GxE':   'PRS\u2082G+GxE',
    'PRS_GxGxE':   'PRS\u2082GxGxE',
}
MODEL_COLORS = {
    'PRS_G':       '#E69F00',   # orange      (main)
    'PRS_G+G':     '#0073A8',   # dark blue   (epi_summed)
    'PRS_GxG':     '#56B4E9',   # sky blue    (epi_product)
    'PRS_G+GxE':   '#006B52',   # dark green  (cardio_summed)
    'PRS_GxGxE':   '#009E73',   # blue-green  (cardio_product)
}

CLUSTER_BG = {
    'HLA autoimmune':        '#fff3f3',
    'Beta-cell (1+2)':       '#f0f7ff',
    'Proinsulin':            '#f0f7ff',
    'Residual glycaemia':    '#f0f7ff',
    'Hyper-ins. secretion':  '#f0f7ff',
    'Liver / lipid metab.':  '#f6fff0',
    'Obesity / body fat':    '#f6fff0',
    'Metabolic syndrome':    '#f6fff0',
    'Lipodystrophy':         '#f6fff0',
}


def load_data(csv_path: str) -> pd.DataFrame:
    raw = pd.read_csv(csv_path)

    # Keep only mapped cohorts and cluster keys
    raw = raw[raw['cohort'].isin(COHORT_MODEL)]
    raw = raw[raw['cluster_key'].isin(CLUSTER_MERGE)]

    raw = raw.copy()
    raw['model']   = raw['cohort'].map(COHORT_MODEL)
    raw['cluster'] = raw['cluster_key'].map(CLUSTER_MERGE)

    # Sum weighted_n across merged clusters (cohort_total_w is constant per cohort)
    grouped = (
        raw.groupby(['model', 'cluster', 'subtype'])
           .agg(weighted_n=('weighted_n', 'sum'),
                cohort_total_w=('cohort_total_w', 'first'))
           .reset_index()
    )
    grouped['weight'] = (grouped['weighted_n'] / grouped['cohort_total_w'] * 100).round(1)
    return grouped[grouped['weight'] > 0].copy()


def make_figure(df: pd.DataFrame, out_dir: str):
    n_clusters = len(CLUSTERS)
    n_subtypes = len(SUBTYPES)

    cluster_idx = {c: i for i, c in enumerate(CLUSTERS)}
    subtype_idx = {s: i for i, s in enumerate(SUBTYPES)}
    model_idx   = {m: i for i, m in enumerate(MODELS)}

    offsets = np.linspace(-0.28, 0.28, len(MODELS))
    SCALE   = 9   # marker area = weight * SCALE

    fig, ax = plt.subplots(figsize=(13, 7))

    # Row background shading
    for ci, cname in enumerate(CLUSTERS):
        y = n_clusters - 1 - ci
        ax.axhspan(y - 0.5, y + 0.5,
                   color=CLUSTER_BG[cname], zorder=0, linewidth=0)

    # Light grid
    for i in range(n_clusters):
        ax.axhline(n_clusters - 1 - i, color='#cccccc', lw=0.4, zorder=1)
    for i in range(n_subtypes):
        ax.axvline(i, color='#cccccc', lw=0.4, zorder=1)

    # Bubbles
    for _, row in df.iterrows():
        if row['cluster'] not in cluster_idx or row['subtype'] not in subtype_idx:
            continue
        ci = cluster_idx[row['cluster']]
        si = subtype_idx[row['subtype']]
        mi = model_idx[row['model']]

        x = si + offsets[mi]
        y = n_clusters - 1 - ci

        ax.scatter(x, y,
                   s=row['weight'] * SCALE,
                   color=MODEL_COLORS[row['model']],
                   alpha=0.82,
                   edgecolors='white',
                   linewidths=0.6,
                   zorder=3)

    # Axes
    ax.set_xticks(range(n_subtypes))
    ax.set_xticklabels(SUBTYPES, fontsize=12, fontweight='bold')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(list(reversed(CLUSTERS)), fontsize=10)

    ax.set_xlim(-0.6, n_subtypes - 0.4)
    ax.set_ylim(-0.6, n_clusters - 0.4)
    ax.tick_params(bottom=False, top=False, left=False)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Biological pathway separators
    immune_end    = n_clusters - 1 - 0   # HLA autoimmune row
    betacell_end  = n_clusters - 1 - 4   # Hyper-ins. secretion row

    ax.axhline(immune_end - 0.5,   color='#aaaaaa', lw=1.2, ls='--', zorder=2)
    ax.axhline(betacell_end - 0.5, color='#aaaaaa', lw=1.2, ls='--', zorder=2)

    ax.text(-0.62, immune_end,                            'Immune',
            va='center', ha='right', fontsize=8.5, color='#c0392b', style='italic')
    ax.text(-0.62, (immune_end - 1 + betacell_end) / 2,  '\u03b2-cell / glycaemia',
            va='center', ha='right', fontsize=8.5, color='#2980b9', style='italic')
    ax.text(-0.62, (betacell_end - 1) / 2,               'Metabolic',
            va='center', ha='right', fontsize=8.5, color='#27ae60', style='italic')

    # Legends
    model_handles = [
        mpatches.Patch(color=MODEL_COLORS[m], label=MODEL_LABELS[m]) for m in MODELS
    ]
    leg1 = fig.legend(handles=model_handles,
                      title='Model',
                      title_fontsize=9,
                      loc='upper right',
                      bbox_to_anchor=(0.99, 0.97),
                      frameon=True,
                      fontsize=8.5,
                      edgecolor='#cccccc',
                      handlelength=1.2)

    size_handles = [
        mlines.Line2D([], [], marker='o', color='gray', alpha=0.75,
                      markersize=np.sqrt(w * SCALE),
                      linestyle='None', label=f'{w}%')
        for w in [5, 15, 30, 50]
    ]
    fig.legend(handles=size_handles,
               title='Weight (%)',
               title_fontsize=9,
               loc='lower right',
               bbox_to_anchor=(0.99, 0.08),
               frameon=True,
               fontsize=9,
               edgecolor='#cccccc',
               handletextpad=1.4)

    ax.set_title(
        'Functional cluster \u2192 T2D subtype: weighted contributions by PRS model',
        fontsize=11, pad=18, loc='left'
    )
    fig.text(0.01, -0.02,
             'Bubble area \u221d weighted % of model features attributed to that cluster\u2013subtype pair.\n'
             'PRS\u2082G+GxE distributes evenly across all subtypes; '
             'PRS\u2082GxGxE is maximally polarised (HLA\u2192SAID, lipodystrophy\u2192SIRD).',
             fontsize=7.5, color='#555555', va='top')

    fig.subplots_adjust(left=0.18, right=0.65, top=0.88, bottom=0.08)

    os.makedirs(out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        out = os.path.join(out_dir, f'cluster_subtype_dotmatrix.{ext}')
        fig.savefig(out, dpi=180, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--data',
        default='results/type2Diabetes/subtypeResearch/enrichment/cluster_subtype_cohort_weighted.csv',
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
