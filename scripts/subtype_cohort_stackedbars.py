#!/usr/bin/env python3
"""
subtype_cohort_stackedbars.py
=============================
Replaces the pie-chart panel (Figure 4B) with two complementary stacked-bar
views that directly translate from the dot-matrix panel (Figure 4A).

Panel B1 — "Subtype totals" (1 × 5 bar chart)
  One bar per T2D subtype.  Bar height = feature count (a meaningful absolute
  quantity, unlike pie percentages that always sum to 100 %).  Each bar is
  colour-stacked by PRS model so the partition across models is visible; the
  hatched portion of each segment is the share of features whose direction
  opposes the cluster's expected direction.

Panel B2 — "Cluster × subtype grid" (cluster rows × subtype columns)
  A small horizontal stacked bar per (cluster, subtype) cell — the same 2-D
  layout as the dot matrix.  Bar length = feature count for that cell; colour
  stack = PRS model; hatching = discordant share.  Empty cells are blank.
  This directly
  shows *which model* drives each dot in the dot-matrix and makes cross-
  subtype counts comparable (unlike pies).

Usage
-----
python scripts/subtype_cohort_stackedbars.py \\
    [--assignments PATH]   # default: results/.../subtype_feature_assignments.csv
    [--out_dir     DIR]    # default: results/.../enrichment/plots
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# ── Identity maps (shared with dotmatrix / pies scripts) ─────────────────────
COHORT_MODEL = {
    'main':           'PRS_G',
    'epi_summed':     'PRS_G+G',
    'epi_product':    'PRS_GxG',
    'cardio_summed':  'PRS_G+GxE',
    'cardio_product': 'PRS_GxGxE',
}

MODEL_LABELS = {
    'PRS_G':     'PRS₂G  (additive)',
    'PRS_G+G':   'PRS₂G+G (G burden)',
    'PRS_GxG':   'PRS₂GxG (epistasis)',
    'PRS_G+GxE': 'PRS₂G+GxE',
    'PRS_GxGxE': 'PRS₂GxGxE',
}

MODEL_COLORS = {
    'PRS_G':     '#E69F00',
    'PRS_G+G':   '#0073A8',
    'PRS_GxG':   '#56B4E9',
    'PRS_G+GxE': '#006B52',
    'PRS_GxGxE': '#009E73',
}

MODELS = ['PRS_G', 'PRS_G+G', 'PRS_GxG', 'PRS_G+GxE', 'PRS_GxGxE']

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

SUBTYPE_FULL = {
    'SAID': 'Severe Autoimmune\nDiabetes',
    'SIDD': 'Severe Insulin-\nDeficient Diabetes',
    'SIRD': 'Severe Insulin-\nResistant Diabetes',
    'MOD':  'Mild Obesity-\nRelated Diabetes',
    'MARD': 'Mild Age-Related\nDiabetes',
}

# Concordance is a DIRECTION (+1 agrees with the cluster's expected direction,
# -1 opposes it, 0 no direction available), not a confidence.  Stacked bars
# need a non-negative length, so bar length is the feature COUNT and the
# discordant share is drawn as a hatched overlay.  An earlier version mapped
# {+1: 1.0, 0: 0.5, -1: 0.0}, which silently dropped discordant features and
# so hid exactly the cases concordance exists to record.

# Set from the data in main(): the hatching note is only shown when at least
# one feature actually runs against its cluster's expected direction, so the
# titles do not advertise an encoding the figure never uses.
DISC_NOTE = ''


# Background shading per cluster (matches dotmatrix)
CLUSTER_BG = {
    'HLA autoimmune':       '#fde8e8',
    'Beta-cell (1+2)':      '#e4eff9',
    'Proinsulin':           '#e4eff9',
    'Residual glycaemia':   '#e4eff9',
    'Hyper-ins. secretion': '#e4eff9',
    'Liver / lipid metab.': '#e8f5e3',
    'Obesity / body fat':   '#e8f5e3',
    'Metabolic syndrome':   '#e8f5e3',
    'Lipodystrophy':        '#e8f5e3',
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(csv_path: str):
    """
    Returns two DataFrames:
      subtype_df  — (subtype, model, weighted_n, n_discordant)   for Panel B1
      cell_df     — (subtype, cluster, model, weighted_n, n_discordant)  for B2

    weighted_n is the feature count (bar length); n_discordant is how many of
    those features run against the cluster's expected direction, drawn as a
    hatched overlay so the disagreement is visible rather than discarded.
    """
    raw = pd.read_csv(csv_path)
    raw = raw[raw['cohort'].isin(COHORT_MODEL)].copy()

    # Drop direct-evidence rows
    raw = raw[~raw['cluster_key'].fillna('').str.startswith('direct__')].copy()

    if 'concordance' in raw.columns:
        raw['_disc'] = (raw['concordance'].fillna(0.0) < 0).astype(int)
    else:
        raw['_disc'] = 0

    dedup = raw.drop_duplicates(subset=['feature', 'cohort', 'cluster_key', 'subtype'])
    dedup = dedup.copy()
    dedup['model'] = dedup['cohort'].map(COHORT_MODEL)

    # Panel B2: per (subtype, cluster, model) — filtered to displayed clusters
    dedup2 = dedup[dedup['cluster_key'].isin(CLUSTER_MERGE)].copy()
    dedup2['cluster'] = dedup2['cluster_key'].map(CLUSTER_MERGE)
    cell_df = (
        dedup2.groupby(['subtype', 'cluster', 'model'])
              .agg(weighted_n=('feature', 'nunique'),
                   n_discordant=('_disc', 'sum'))
              .reset_index()
    )
    cell_df = cell_df[cell_df['weighted_n'] > 0].copy()

    # Panel B1: totals derived from the same filtered clusters so column totals
    # in B1 equal the sum of cell values in B2's column.
    subtype_df = (
        cell_df.groupby(['subtype', 'model'])[['weighted_n', 'n_discordant']]
               .sum()
               .reset_index()
    )

    return subtype_df, cell_df


# ── Panel B1: Subtype-level stacked bar chart ─────────────────────────────────

def draw_subtype_bars(subtype_df: pd.DataFrame, out_path: str):
    """
    One vertical stacked bar per subtype.  Height = absolute concordance-
    weighted feature count.  Colour = PRS model.
    """
    fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')

    x = np.arange(len(SUBTYPES))
    bar_width = 0.55
    bottoms = np.zeros(len(SUBTYPES))

    for model in MODELS:
        vals = []
        for st in SUBTYPES:
            row = subtype_df[(subtype_df['subtype'] == st) &
                             (subtype_df['model'] == model)]
            vals.append(row['weighted_n'].values[0] if len(row) else 0.0)
        vals = np.array(vals)
        discs = np.array([
            subtype_df.loc[(subtype_df['subtype'] == st) &
                           (subtype_df['model'] == model), 'n_discordant'].sum()
            for st in SUBTYPES
        ], dtype=float)
        ax.bar(x, vals, bar_width, bottom=bottoms,
               color=MODEL_COLORS[model], label=MODEL_LABELS[model],
               edgecolor='white', linewidth=0.6)
        # Discordant share, hatched at the top of the model's segment
        ax.bar(x, discs, bar_width, bottom=bottoms + vals - discs,
               color='none', edgecolor='white', hatch='////',
               linewidth=0.0, zorder=3)
        bottoms += vals

    # Total annotation above each bar
    for xi, total in enumerate(bottoms):
        ax.text(xi, total + 0.5, f'{total:.0f}',
                ha='center', va='bottom', fontsize=8.5, color='#333333')

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f'{s}\n{SUBTYPE_FULL[s]}' for s in SUBTYPES],
        fontsize=8.5,
    )
    ax.set_ylabel('Feature count', fontsize=10)
    ax.set_title(
        f'Features per subtype, partitioned by PRS model{DISC_NOTE}',
        fontsize=11, fontweight='bold', pad=10,
    )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(bottom=False)
    ax.yaxis.grid(True, color='#e0e0e0', zorder=0)
    ax.set_axisbelow(True)

    handles = [mpatches.Patch(color=MODEL_COLORS[m], label=MODEL_LABELS[m])
               for m in MODELS]
    ax.legend(handles=handles, title='PRS Model', title_fontsize=9,
              fontsize=8.5, framealpha=0.95, edgecolor='#cccccc',
              loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)

    fig.subplots_adjust(right=0.72)
    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved B1 → {out_path}')


# ── Panel B2: Cluster × subtype grid — single-axes, matches dot-matrix layout ──

# Width of one subtype column in data units (bar + gap)
COL_WIDTH = 1.0
# Fraction of COL_WIDTH used by the bar; rest is inter-column gap
BAR_FILL  = 0.72

# Section row-index boundaries (same as dot-matrix)
# immune=0, beta-cell=1-4, metabolic=5-8
SECTION_SEP_ROWS = [0.5, 4.5]   # dashed lines drawn at these y values

SECTION_LABELS = [
    (0,   0,   'Immune',         '#c0392b'),
    (1,   4,   'β-cell /\nglycaemia', '#2477b5'),
    (5,   8,   'Metabolic',      '#2a7a2a'),
]


def _draw_grid_ax(ax, cell_df: pd.DataFrame):
    """
    Populate a single axes with the cluster × subtype stacked-bar grid.
    Rows = clusters (top→bottom), columns = subtypes.
    Returns global_max used for the x-scale (per column).
    """
    n_c = len(CLUSTERS)
    n_s = len(SUBTYPES)

    cell_totals = cell_df.groupby(['subtype', 'cluster'])['weighted_n'].sum()
    global_max = cell_totals.max() if len(cell_totals) else 1.0

    bar_w  = COL_WIDTH * BAR_FILL          # bar height in y-units
    bar_h  = 0.70                          # bar thickness in y-units

    # ── background row shading ────────────────────────────────────────────────
    for ci, cluster in enumerate(CLUSTERS):
        y_center = n_c - 1 - ci
        color = CLUSTER_BG.get(cluster, '#f9f9f9')
        ax.axhspan(y_center - 0.5, y_center + 0.5, color=color,
                   zorder=0, linewidth=0)

    # ── white row-anchor lines drawn above background, below bars ────────────
    for ci in range(n_c):
        ax.axhline(n_c - 1 - ci, color='white', lw=0.8, zorder=2)

    # ── stacked bars ──────────────────────────────────────────────────────────
    for ci, cluster in enumerate(CLUSTERS):
        y_center = n_c - 1 - ci
        for si, subtype in enumerate(SUBTYPES):
            sub = cell_df[
                (cell_df['cluster'] == cluster) &
                (cell_df['subtype'] == subtype)
            ]
            if sub.empty:
                continue

            # x origin for this column's bar, centred in the column cell
            x_left = si * COL_WIDTH + (COL_WIDTH - bar_w) / 2
            cursor = x_left
            total  = 0.0
            for model in MODELS:
                row = sub[sub['model'] == model]
                val = row['weighted_n'].values[0] if len(row) else 0.0
                if val <= 0:
                    continue
                disc = row['n_discordant'].values[0] if len(row) else 0.0
                seg_w = val / global_max * bar_w
                rect = Rectangle(
                    (cursor, y_center - bar_h / 2),
                    seg_w, bar_h,
                    facecolor=MODEL_COLORS[model],
                    edgecolor='none',
                    linewidth=0,
                    antialiased=False,
                    zorder=3,
                )
                ax.add_patch(rect)
                if disc > 0:
                    disc_w = disc / global_max * bar_w
                    ax.add_patch(Rectangle(
                        (cursor + seg_w - disc_w, y_center - bar_h / 2),
                        disc_w, bar_h,
                        facecolor='none', edgecolor='white', hatch='////',
                        linewidth=0, zorder=4,
                    ))
                cursor += seg_w
                total  += val

            # Total label to the right of the bar
            ax.text(si * COL_WIDTH + COL_WIDTH - 0.03, y_center,
                    f'{total:.0f}',
                    ha='right', va='center', fontsize=7, color='#333333',
                    zorder=4)

    # ── white row-anchor lines redrawn on top of bars ────────────────────────
    for ci in range(n_c):
        ax.axhline(n_c - 1 - ci, color='white', lw=0.8, zorder=5)

    # ── dashed vertical separators between subtype columns ────────────────────
    for si in range(1, n_s):
        ax.axvline(si * COL_WIDTH, color='#aaaaaa', lw=0.8,
                   ls='--', zorder=2)

    # ── dashed horizontal lines at major section boundaries only ─────────────
    # immune (row 0) | β-cell/glycaemia (rows 1-4) | metabolic (rows 5-8)
    for ci_boundary in [1, 5]:   # first row of each new section
        y_sep = n_c - 1 - ci_boundary + 0.5
        ax.axhline(y_sep, color='#aaaaaa', lw=1.0, ls='--', zorder=2)

    # ── axes limits and ticks ─────────────────────────────────────────────────
    ax.set_xlim(0, n_s * COL_WIDTH)
    ax.set_ylim(-0.55, n_c - 0.45)

    # Column (subtype) tick marks centred in each column cell
    ax.set_xticks([si * COL_WIDTH + COL_WIDTH / 2 for si in range(n_s)])
    ax.set_xticklabels(SUBTYPES, fontsize=11, fontweight='bold')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(bottom=False, top=False, left=False)

    # Row (cluster) labels on y-axis
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(list(reversed(CLUSTERS)), fontsize=9.5)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # ── section italic labels in left margin (further out to avoid row labels) ─
    for r_start, r_end, label, color in SECTION_LABELS:
        y_top = n_c - 1 - r_start
        y_bot = n_c - 1 - r_end
        y_mid = (y_top + y_bot) / 2
        ax.text(-3.2, y_mid, label,
                ha='center', va='center', fontsize=8.5,
                color=color, style='italic', clip_on=False,
                rotation=90)

    return global_max




def draw_grid_bars(cell_df: pd.DataFrame, out_path: str):
    """Single-axes dot-matrix-style grid with stacked bars per cell."""
    fig, ax = plt.subplots(figsize=(13.5, 8.5), facecolor='white')

    _draw_grid_ax(ax, cell_df)

    handles = [mpatches.Patch(color=MODEL_COLORS[m], label=MODEL_LABELS[m])
               for m in MODELS]
    fig.legend(handles=handles, title='PRS Model', title_fontsize=9,
               fontsize=8.5, framealpha=0.97, edgecolor='#cccccc',
               loc='upper right', bbox_to_anchor=(1.00, 0.97))

    ax.set_title(
        f'Feature count by cluster × subtype, by PRS model{DISC_NOTE}',
        fontsize=11, fontweight='bold', pad=18,
    )

    fig.subplots_adjust(left=0.18, right=0.78)
    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved B2 → {out_path}')


# ── Combined two-panel figure ─────────────────────────────────────────────────

def draw_combined(subtype_df: pd.DataFrame, cell_df: pd.DataFrame,
                  out_path: str):
    """
    Combined B1 (top) + B2 (bottom) in one figure for a direct Figure 4B
    replacement.
    """
    fig = plt.figure(figsize=(13, 13), facecolor='white')
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 2.2], hspace=0.35)

    # ── B1: subtype totals ────────────────────────────────────────────────────
    ax_top = fig.add_subplot(gs[0])

    x = np.arange(len(SUBTYPES))
    bar_width = 0.55
    bottoms = np.zeros(len(SUBTYPES))

    for model in MODELS:
        vals = np.array([
            subtype_df.loc[(subtype_df['subtype'] == st) &
                           (subtype_df['model'] == model), 'weighted_n'].sum()
            for st in SUBTYPES
        ])
        discs = np.array([
            subtype_df.loc[(subtype_df['subtype'] == st) &
                           (subtype_df['model'] == model), 'n_discordant'].sum()
            for st in SUBTYPES
        ], dtype=float)
        ax_top.bar(x, vals, bar_width, bottom=bottoms,
                   color=MODEL_COLORS[model], label=MODEL_LABELS[model],
                   edgecolor='white', linewidth=0.6)
        ax_top.bar(x, discs, bar_width, bottom=bottoms + vals - discs,
                   color='none', edgecolor='white', hatch='////',
                   linewidth=0.0, zorder=3)
        bottoms += vals

    for xi, total in enumerate(bottoms):
        ax_top.text(xi, total + 0.3, f'{total:.0f}',
                    ha='center', va='bottom', fontsize=8.5, color='#333333')

    ax_top.set_xticks(x)
    ax_top.set_xticklabels(
        [f'{s}\n{SUBTYPE_FULL[s]}' for s in SUBTYPES], fontsize=8.5
    )
    ax_top.set_ylabel('Feature count', fontsize=9.5)
    ax_top.set_title(f'B1  Features per subtype, partitioned by PRS model{DISC_NOTE}',
                     fontsize=10.5, fontweight='bold', loc='left', pad=6)
    for spine in ['top', 'right']:
        ax_top.spines[spine].set_visible(False)
    ax_top.tick_params(bottom=False)
    ax_top.yaxis.grid(True, color='#e0e0e0', zorder=0)
    ax_top.set_axisbelow(True)

    handles = [mpatches.Patch(color=MODEL_COLORS[m], label=MODEL_LABELS[m])
               for m in MODELS]
    ax_top.legend(handles=handles, title='PRS Model', title_fontsize=8.5,
                  fontsize=8, framealpha=0.95, edgecolor='#cccccc',
                  loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)

    # ── B2: single-axes dot-matrix-style grid ────────────────────────────────
    ax_grid = fig.add_subplot(gs[1])
    _draw_grid_ax(ax_grid, cell_df)
    ax_grid.set_title(
        f'B2  Per-cluster features (rows) × subtype (columns), by PRS model{DISC_NOTE}',
        fontsize=10.5, fontweight='bold', loc='left', pad=18,
    )

    fig.subplots_adjust(left=0.18, right=0.78)
    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved combined → {out_path}')


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        '--assignments',
        default='results/type2Diabetes/combinedAnalysis/scores/'
                'functionalEnrichment/subtype_feature_assignments.csv',
    )
    p.add_argument(
        '--out_dir',
        default='results/type2Diabetes/subtypeResearch/enrichment/plots',
    )
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f'Loading {args.assignments}')
    subtype_df, cell_df = load_data(args.assignments)
    global DISC_NOTE
    if cell_df['n_discordant'].sum() > 0:
        DISC_NOTE = ' (hatched = discordant)'

    draw_subtype_bars(
        subtype_df,
        os.path.join(args.out_dir, 'subtype_cohort_bars_B1.png'),
    )
    draw_grid_bars(
        cell_df,
        os.path.join(args.out_dir, 'subtype_cohort_bars_B2_grid.png'),
    )
    draw_combined(
        subtype_df, cell_df,
        os.path.join(args.out_dir, 'subtype_cohort_bars_combined.png'),
    )
    print('\nDone.')


if __name__ == '__main__':
    main()
