#!/usr/bin/env python3
"""
subtype_cohort_pies.py
======================
One donut-style pie chart per T2D subtype showing the relative concordance-
weighted cohort contribution to each subtype's feature pool.

Concordance weighting (per Udler 2019 / Suzuki 2024 cluster phenotype profiles):
  +1  concordant  → full weight  (1.0)
   0  neutral     → half weight  (0.5)   genomic loci without a clinical direction map
  -1  discordant  → zero weight  (0.0)   feature direction opposes cluster expectation

Two output files are written per subtype:
  subtype_pie_{SUBTYPE}_labelled.png   – wedge %, cohort name, raw n
  subtype_pie_{SUBTYPE}_clean.png      – no text labels (for manual Inkscape annotation)

Usage
-----
python subtype_cohort_pies.py \\
    --assignments results/.../subtype_feature_assignments.csv \\
    --out_dir     results/.../enrichment/plots
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── palette ──────────────────────────────────────────────────────────────────
COHORT_COLORS = {
    'epi_product':    '#56B4E9',   # sky blue   (Okabe-Ito)
    'epi_summed':     '#0073A8',   # darker blue
    'cardio_product': '#009E73',   # bluish green (Okabe-Ito)
    'cardio_summed':  '#006B52',   # darker green
    'main':           '#E69F00',   # orange (Okabe-Ito)
}

COHORT_LABELS = {
    'cardio_product': 'Cardio Product',
    'cardio_summed':  'Cardio Summed',
    'epi_product':    'Epi Product',
    'epi_summed':     'Epi Summed',
    'main':           'Main',
}

SUBTYPE_ORDER = ['SAID', 'SIDD', 'SIRD', 'MOD', 'MARD']

SUBTYPE_FULL = {
    'SAID': 'Severe Autoimmune\nDiabetes (SAID)',
    'SIDD': 'Severe Insulin-Deficient\nDiabetes (SIDD)',
    'SIRD': 'Severe Insulin-Resistant\nDiabetes (SIRD)',
    'MOD':  'Mild Obesity-Related\nDiabetes (MOD)',
    'MARD': 'Mild Age-Related\nDiabetes (MARD)',
}

# ── concordance weight mapping ────────────────────────────────────────────────
CONC_WEIGHT = {1.0: 1.0, 0.0: 0.5, -1.0: 0.0}


def compute_weighted_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute concordance-weighted feature count per (subtype, cohort).

    Each unique (feature, cohort, cluster_key, subtype) row contributes
    weight = CONC_WEIGHT[concordance] to its (subtype, cohort) bucket.
    De-duplication is done first so a feature appearing in two clusters for
    the same subtype is counted once per cluster (as in the assignment table).
    """
    has_conc = 'concordance' in df.columns
    if has_conc:
        df = df.copy()
        df['weight'] = df['concordance'].map(CONC_WEIGHT).fillna(0.5)
    else:
        df = df.copy()
        df['weight'] = 0.5  # treat all as neutral if column absent

    # Sum weights per subtype × cohort
    agg = (df.groupby(['subtype', 'cohort'])['weight']
             .sum()
             .reset_index(name='weighted_n'))

    # Also keep raw count for label annotation
    raw = (df.groupby(['subtype', 'cohort'])
             .size()
             .reset_index(name='raw_n'))

    return agg.merge(raw, on=['subtype', 'cohort'])


def draw_pie(subtype: str, group: pd.DataFrame, out_path: str,
             show_labels: bool = True):
    """Draw a single donut-pie for one subtype and save to out_path."""
    cohort_order = [c for c in list(COHORT_COLORS) if c in group['cohort'].values]
    group = group.set_index('cohort').reindex(cohort_order).dropna(subset=['weighted_n'])
    group = group[group['weighted_n'] > 0]

    if group.empty:
        print(f'  [{subtype}] no positive-weight cohorts – skipping')
        return

    sizes  = group['weighted_n'].values
    raws   = group['raw_n'].values.astype(int)
    colors = [COHORT_COLORS[c] for c in group.index]
    labels = [COHORT_LABELS.get(c, c) for c in group.index]

    fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    ax.set_aspect('equal')

    total = sizes.sum()

    # ── donut wedges ──────────────────────────────────────────────────────────
    wedges, _ = ax.pie(
        sizes,
        colors=colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=0.52, edgecolor='white', linewidth=2.0),
    )

    # ── centre annotation ─────────────────────────────────────────────────────
    acronym_y = 0.08 if show_labels else 0.0
    ax.text(0, acronym_y, subtype,
            ha='center', va='center',
            fontsize=20, fontweight='bold', color='#333333')
    if show_labels:
        ax.text(0, -0.18, SUBTYPE_FULL.get(subtype, subtype).replace('\n', ' '),
                ha='center', va='center',
                fontsize=7.5, color='#777777', style='italic')

    # ── wedge labels (labelled version only) ─────────────────────────────────
    if show_labels:
        for wedge, sz, raw, lbl in zip(wedges, sizes, raws, labels):
            pct = 100.0 * sz / total
            if pct < 3:          # skip tiny slices to avoid clutter
                continue
            # Midpoint angle of wedge (matplotlib measures CCW from +x, but we
            # set startangle=90 CCW=False so we convert back)
            ang = np.radians((wedge.theta1 + wedge.theta2) / 2)
            r_lbl = 0.72          # just outside the donut ring
            x = r_lbl * np.cos(ang)
            y = r_lbl * np.sin(ang)
            txt = f'{lbl}\n{pct:.1f}%\n(n={raw})'
            ax.text(x, y, txt,
                    ha='center', va='center',
                    fontsize=8.0, color='#222222',
                    bbox=dict(boxstyle='round,pad=0.25', fc='white',
                              ec='#cccccc', alpha=0.85, linewidth=0.8))

    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.axis('off')

    # ── legend (always shown) ─────────────────────────────────────────────────
    handles = [mpatches.Patch(color=COHORT_COLORS[c], label=COHORT_LABELS.get(c, c))
               for c in cohort_order if c in group.index]
    ax.legend(handles=handles,
              loc='lower center', bbox_to_anchor=(0.5, -0.07),
              ncol=3, fontsize=8.5, framealpha=0.95, edgecolor='#cccccc',
              title='Cohort', title_fontsize=8.5)

    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [{subtype}]  →  {out_path}')


def draw_panel(weighted: pd.DataFrame, out_path: str, show_labels: bool = True):
    """
    Draw all subtypes as a single 1×5 panel figure and save to out_path.
    This is the primary output; individual files are also saved.
    """
    subtypes_present = [s for s in SUBTYPE_ORDER if s in weighted['subtype'].values]
    n = len(subtypes_present)

    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.8), facecolor='white')
    if n == 1:
        axes = [axes]

    total_per_subtype = weighted.groupby('subtype')['weighted_n'].sum()

    for ax, subtype in zip(axes, subtypes_present):
        group = weighted[weighted['subtype'] == subtype].copy()
        cohort_order = [c for c in list(COHORT_COLORS) if c in group['cohort'].values]
        group = group.set_index('cohort').reindex(cohort_order).dropna(subset=['weighted_n'])
        group = group[group['weighted_n'] > 0]

        if group.empty:
            ax.axis('off')
            continue

        sizes  = group['weighted_n'].values
        raws   = group['raw_n'].values.astype(int)
        colors = [COHORT_COLORS[c] for c in group.index]
        labels = [COHORT_LABELS.get(c, c) for c in group.index]
        total  = sizes.sum()

        wedges, _ = ax.pie(
            sizes,
            colors=colors,
            startangle=90,
            counterclock=False,
            wedgeprops=dict(width=0.52, edgecolor='white', linewidth=2.0),
        )

        # Centre text
        acronym_y = 0.10 if show_labels else 0.0
        ax.text(0, acronym_y, subtype,
                ha='center', va='center',
                fontsize=16, fontweight='bold', color='#333333',
                transform=ax.transData)
        if show_labels:
            ax.text(0, -0.15,
                    SUBTYPE_FULL.get(subtype, subtype),
                    ha='center', va='center',
                    fontsize=6.5, color='#888888', style='italic',
                    transform=ax.transData)

        if show_labels:
            for wedge, sz, raw, lbl in zip(wedges, sizes, raws, labels):
                pct = 100.0 * sz / total
                if pct < 4:
                    continue
                ang = np.radians((wedge.theta1 + wedge.theta2) / 2)
                r_lbl = 0.74
                x = r_lbl * np.cos(ang)
                y = r_lbl * np.sin(ang)
                txt = f'{pct:.0f}%\n(n={raw})'
                ax.text(x, y, txt,
                        ha='center', va='center', fontsize=7.5, color='#222222',
                        bbox=dict(boxstyle='round,pad=0.2', fc='white',
                                  ec='#cccccc', alpha=0.82, linewidth=0.7))

        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.set_aspect('equal')
        ax.axis('off')

    # Shared legend below the panel
    handles = [mpatches.Patch(color=v, label=COHORT_LABELS.get(k, k))
               for k, v in COHORT_COLORS.items()
               if k in weighted['cohort'].values]
    fig.legend(handles=handles,
               loc='lower center', bbox_to_anchor=(0.5, -0.02),
               ncol=len(handles), fontsize=9.5, framealpha=0.97,
               edgecolor='#cccccc', title='PRS Cohort', title_fontsize=9.5)

    tag = '(concordance-weighted feature count)' if show_labels else ''
    fig.suptitle(f'Cohort Enrichment per T2D Subtype  {tag}',
                 fontsize=13, fontweight='bold', y=1.01)

    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Panel  →  {out_path}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',
                   default='results/type2Diabetes/subtypeResearch/enrichment/'
                           'subtype_feature_assignments.csv')
    p.add_argument('--out_dir',
                   default='results/type2Diabetes/subtypeResearch/enrichment/plots')
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f'Loading {args.assignments}')
    df = pd.read_csv(args.assignments)

    # Drop direct-evidence rows (cluster_key starts with 'direct__') so they
    # don't double-count features already counted via the cluster route
    df = df[~df['cluster_key'].fillna('').str.startswith('direct__')].copy()

    # De-duplicate: one row per unique (feature, cohort, cluster_key, subtype)
    df = df.drop_duplicates(subset=['feature', 'cohort', 'cluster_key', 'subtype'])

    weighted = compute_weighted_counts(df)
    print('\nConcordance-weighted counts:')
    print(weighted.sort_values(['subtype', 'cohort']).to_string(index=False))

    # ── panel outputs ──────────────────────────────────────────────────────────
    draw_panel(weighted,
               os.path.join(args.out_dir, 'subtype_cohort_pies_labelled.png'),
               show_labels=True)
    draw_panel(weighted,
               os.path.join(args.out_dir, 'subtype_cohort_pies_clean.png'),
               show_labels=False)

    # ── individual per-subtype outputs ────────────────────────────────────────
    for subtype in SUBTYPE_ORDER:
        grp = weighted[weighted['subtype'] == subtype].copy()
        if grp.empty:
            continue
        draw_pie(subtype, grp,
                 os.path.join(args.out_dir, f'subtype_pie_{subtype}_labelled.png'),
                 show_labels=True)
        draw_pie(subtype, grp,
                 os.path.join(args.out_dir, f'subtype_pie_{subtype}_clean.png'),
                 show_labels=False)

    print('\nDone.')


if __name__ == '__main__':
    main()
