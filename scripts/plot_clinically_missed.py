#!/usr/bin/env python3
# Keri Multerer 2026
# Figure 6: what ePRS_Comp contributes among individuals clinical measures do
# not flag.
#
# WHAT THIS REPLACES
# The previous Figure 6 plotted NRI between a clinical threshold rule and a PRS
# rule. Those are not nested models and flagged very different fractions of the
# cohort -- under 1% for glucose against 64.4% for the composite -- so the
# reclassification counts were set by that mismatch rather than by
# discrimination. It also defined "clinically low risk" with glucose at
# 12 mmol/L, above the diagnostic cut-point, which left 99.4% of participants
# in the low-risk group: there was effectively no subgroup.
#
# WHAT IS PLOTTED
# Both scores are evaluated at the SAME flag rate inside each subgroup, so the
# comparison is at a matched classification budget and neither can look better
# merely by flagging more people.
#
#   Panel A  precision (the share of flagged individuals who are cases) for
#            each score, against that subgroup's own case rate. The gap between
#            bar and reference line is the enrichment over flagging at random.
#   Panel B  cases captured at that budget, which is the same quantity in
#            absolute terms.
#
# Panel A answers "how concentrated is the yield" and Panel B "how many people
# is that", which is the pair a clinician needs: a high precision on a tiny
# subgroup and a large count on a diffuse one are different propositions.
#
# Subgroups are ordered by case rate so the reference line falls monotonically
# and the reader can see precision rising against a falling baseline.

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SCORE_COLOURS = {'ePRS_Comp': '#D55E00', 'PRS_G': '#E69F00'}
SCORE_LABELS = {'ePRS_Comp': '$\\mathrm{ePRS}_{\\mathrm{Comp}}$', 'PRS_G': '$\\mathrm{PRS}_{G}$'}

# Long analysis labels shortened for the axis; the full definitions live in the
# caption and the Methods.
SHORT = {
    'whole holdout': 'Whole cohort',
    'not flagged: Glycated haemoglobin (HbA1c)': 'HbA1c\nnot flagged',
    'not flagged: Body mass index (BMI)': 'BMI\nnot flagged',
    'not flagged: Glucose': 'Glucose\nnot flagged',
    'biomarker-negative': 'No measure\nflagged',
    'model-negative': 'Clinical model\nnot flagged',
}


def draw_single(d, order, tag, out_dir):
    """One axes, precision on the left and case counts on the right.

    The two scales are independent, so their relative heights carry no
    meaning: only each series against its own axis is interpretable. Bars are
    used for precision and markers for counts so the eye does not read one as
    continuous with the other, and the axis labels are colour-matched to the
    elements they govern.
    """
    x = np.arange(len(order))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    ax2 = ax.twinx()

    for i, score in enumerate(['ePRS_Comp', 'PRS_G']):
        sub = d[d.score == score].set_index('label').loc[order]
        off = (i - 0.5) * width
        ax.bar(x + off, sub[f'{tag}_ppv'] * 100, width * 0.92,
               color=SCORE_COLOURS[score], label=f'{SCORE_LABELS[score]} precision',
               edgecolor='white', linewidth=0.6, zorder=3)
        ax2.plot(x + off, sub[f'{tag}_cases_caught'],
                 marker='D' if i == 0 else 'o', markersize=7, linestyle='none',
                 markerfacecolor='white', markeredgewidth=1.8,
                 markeredgecolor=SCORE_COLOURS[score],
                 label=f'{SCORE_LABELS[score]} cases', zorder=5)
        for xi, c in enumerate(sub[f'{tag}_cases_caught']):
            ax2.text(xi + off, c, f'  {int(c)}', ha='left', va='center',
                     fontsize=7.5, color=SCORE_COLOURS[score])

    base = d[d.score == 'ePRS_Comp'].set_index('label').loc[order, 'base_rate'] * 100
    for xi, b in enumerate(base):
        ax.plot([xi - width, xi + width], [b, b], color='#444444', lw=1.4, ls='--',
                zorder=4, label='Case rate in subgroup' if xi == 0 else None)

    ax.set_ylabel('Precision at the top 20 % (%), bars', fontsize=9.5, color='#8a3b00')
    ax2.set_ylabel('Cases captured, markers', fontsize=9.5, color='#333333')
    ax.tick_params(axis='y', colors='#8a3b00')
    ax.set_ylim(0, max(base.max(), (d[f'{tag}_ppv'] * 100).max()) * 1.35)
    ax2.set_ylim(0, d[f'{tag}_cases_caught'].max() * 1.35)

    n = d[d.score == 'ePRS_Comp'].set_index('label').loc[order]
    ax.set_xticks(x)
    ax.set_xticklabels([f'{lab}\nn={int(r.n):,}; {int(r.n_cases):,} cases'
                        for lab, r in zip(order, n.itertuples())], fontsize=8)
    ax.set_title('$\\mathrm{ePRS}_{\\mathrm{Comp}}$ among individuals clinical measures do not flag',
                 fontsize=11, fontweight='bold', loc='left', pad=8)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8.5, frameon=True, edgecolor='#cccccc',
              loc='upper right', ncol=2)
    for sp in ('top',):
        ax.spines[sp].set_visible(False); ax2.spines[sp].set_visible(False)
    ax.yaxis.grid(True, color='#e8e8e8', zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(bottom=False)

    fig.text(0.01, -0.06,
             'Both scores flag the top 20 % of each subgroup, so they are compared at the same classification '
             'budget. Bars read against the left axis and markers against the right;\nthe two scales are '
             'independent, so the height of a bar relative to a marker has no meaning. Subgroups are ordered '
             'by case rate.',
             fontsize=7.5, color='#555555', va='top')
    fig.subplots_adjust(left=0.08, right=0.92, top=0.9, bottom=0.2)
    for ext in ('png', 'pdf'):
        out = os.path.join(out_dir, f'clinically_missed_subgroups_single.{ext}')
        fig.savefig(out, dpi=180, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)


SUBGROUPS_COMPACT = [
    ('whole holdout',                              'Whole cohort'),
    ('not flagged: Glycated haemoglobin (HbA1c)', 'HbA1c\nnot flagged'),
    ('not flagged: Glucose',                       'Glucose\nnot flagged'),
    ('not flagged: Body mass index (BMI)',          'BMI\nnot flagged'),
    ('biomarker-negative',                          'No measure\nflagged'),
]


def draw_compact(decile_csv, out_dir):
    """Precision by score decile for 5 subgroups in one figure.

    Layout per subgroup: ePRS_Comp deciles 10->1 (touching), small gap,
    PRS_G deciles 10->1 (touching), larger gap to next subgroup.
    A light dashed line marks each subgroup's overall case rate.
    Bars above that line are annotated with additional cases above baseline.
    """
    dec = pd.read_csv(decile_csv)

    BLOCK = 10      # deciles per score
    BW = 0.46       # bar width; adjacent bars within a block touch (center spacing = BW)
    SGAP = 0.10     # small gap between ePRS_Comp block and PRS_G block
    GGAP = 0.28     # gap between subgroups — slightly wider than the separator line
    SPAN = 2 * BLOCK * BW + SGAP + GGAP   # x-units per subgroup

    def block_start(gi, score_idx):
        """Left edge of the score block (0=ePRS_Comp, 1=PRS_G) in group gi."""
        base = gi * SPAN
        return base if score_idx == 0 else base + BLOCK * BW + SGAP

    def bar_center(gi, score_idx, d):
        """Centre x of the bar for decile d (1-10), descending: d=10 -> pos 0."""
        return block_start(gi, score_idx) + (BLOCK - d) * BW + BW / 2

    fig, ax = plt.subplots(figsize=(14, 5.0))
    all_xticks, all_xticklabels = [], []

    for gi, (sub_key, sub_label) in enumerate(SUBGROUPS_COMPACT):
        gdata = dec[dec.subgroup == sub_key].copy()
        if gdata.empty:
            continue

        comp_d = gdata[gdata.score == 'ePRS_Comp']
        base = comp_d['n_cases'].sum() / comp_d['n'].sum()
        total_n = int(comp_d['n'].sum())
        total_cases = int(comp_d['n_cases'].sum())

        for si, score in enumerate(['ePRS_Comp', 'PRS_G']):
            sd = gdata[gdata.score == score].sort_values('decile', ascending=False)
            xs = np.array([bar_center(gi, si, d) for d in sd['decile']])
            ax.bar(xs, sd['case_rate'] * 100, BW,
                   color=SCORE_COLOURS[score],
                   label=SCORE_LABELS[score] if gi == 0 else '_nolegend_',
                   edgecolor='none', linewidth=0, zorder=3)
            # Additional cases above subgroup baseline
            for xi, row in zip(xs, sd.itertuples()):
                excess = row.n * (row.case_rate - base)
                if excess > 0.5:
                    ax.annotate(f'+{int(round(excess))}',
                                (xi, row.case_rate * 100),
                                textcoords='offset points', xytext=(0, 2),
                                ha='center', va='bottom', fontsize=5.5,
                                color='#333333', annotation_clip=False)

            # Decile tick labels below each block
            for d in range(10, 0, -1):
                all_xticks.append(bar_center(gi, si, d))
                all_xticklabels.append(str(d))

            # Score label centred below the block
            block_mid = block_start(gi, si) + BLOCK * BW / 2
            ax.annotate(SCORE_LABELS[score],
                        xy=(block_mid, 0), xycoords='data',
                        xytext=(0, -26), textcoords='offset points',
                        ha='center', va='top', fontsize=7.5,
                        color=SCORE_COLOURS[score], annotation_clip=False,
                        fontweight='bold')

        # Dashed base-rate line spanning both score blocks
        x_lo = block_start(gi, 0)
        x_hi = block_start(gi, 1) + BLOCK * BW
        ax.hlines(base * 100, x_lo, x_hi, colors='#888888', lw=0.9, ls='--',
                  zorder=4, label='Subgroup case rate' if gi == 0 else '_nolegend_')

        # Subgroup label + n centred below both blocks
        sg_mid = block_start(gi, 0) + BLOCK * BW + SGAP / 2
        ax.annotate(f'{sub_label}\nn={total_n:,}; {total_cases:,} cases',
                    xy=(sg_mid, 0), xycoords='data',
                    xytext=(0, -48), textcoords='offset points',
                    ha='center', va='top', fontsize=8, color='#333333',
                    annotation_clip=False)

        # Vertical separator between subgroups
        if gi < len(SUBGROUPS_COMPACT) - 1:
            sep = block_start(gi, 1) + BLOCK * BW + GGAP / 2
            ax.axvline(sep, color='#cccccc', lw=1.0, zorder=2)

    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xticklabels, fontsize=6)
    ax.set_ylabel('Precision (% of decile who are cases)', fontsize=9.5)
    ax.set_ylim(0, dec['case_rate'].max() * 100 * 1.32)
    n_groups = len(SUBGROUPS_COMPACT)
    ax.set_xlim(-BW * 0.3,
                block_start(n_groups - 1, 1) + BLOCK * BW + BW * 0.3)
    ax.set_title(
        '$\\mathrm{ePRS}_{\\mathrm{Comp}}$ vs $\\mathrm{PRS}_{G}$'
        '  —  precision by score decile in each subgroup'
        '  (decile 10 = highest score)',
        fontsize=10.5, fontweight='bold', loc='left', pad=8)
    ax.legend(fontsize=8.5, frameon=True, edgecolor='#cccccc', loc='upper right')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.yaxis.grid(True, color='#e8e8e8', zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(bottom=False)

    fig.text(0.01, -0.03,
             'Precision is the fraction of individuals in that score decile who are cases. '
             'Decile boundaries are computed separately for each score within each subgroup. '
             'The dashed line is the subgroup\'s overall case rate; '
             'bars above it are annotated with the number of excess cases\n'
             'relative to random selection from that subgroup. '
             'Decile 10 = highest-scoring individuals.',
             fontsize=7.5, color='#555555', va='top')
    fig.subplots_adjust(left=0.04, right=0.99, top=0.90, bottom=0.38)
    for ext in ('png', 'pdf'):
        out = os.path.join(out_dir, f'clinically_missed_subgroups_compact.{ext}')
        fig.savefig(out, dpi=180, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results_csv', required=True,
                    help='eprs_in_clinically_missed.csv')
    ap.add_argument('--decile_csv', default=None,
                    help='eprs_in_clinically_missed_deciles.csv '
                         '(required for --layout compact)')
    ap.add_argument('--subgroup', default='biomarker-negative',
                    help='subgroup key to show in the compact decile plot')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--flag_rate', type=float, default=0.20)
    ap.add_argument('--layout', choices=['two_panel', 'single', 'compact'],
                    default='two_panel',
                    help="'single' overlays both quantities on one axes with "
                         "precision on the left and case counts on the right; "
                         "'compact' shows deciles 1-10 with thin bars.")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tag = f'top{int(args.flag_rate * 100)}'
    d = pd.read_csv(args.results_csv)
    d = d[d.subgroup.isin(SHORT)].copy()
    d['label'] = d.subgroup.map(SHORT)

    # Order by case rate: the reference line then falls left to right, so a
    # rising bar against it reads directly as increasing enrichment.
    order = (d[d.score == 'ePRS_Comp']
             .sort_values('base_rate', ascending=False)['label'].tolist())
    if args.layout == 'single':
        draw_single(d, order, tag, args.out_dir)
        return
    if args.layout == 'compact':
        if not args.decile_csv:
            ap.error('--decile_csv is required for --layout compact')
        draw_compact(args.decile_csv, args.out_dir)
        return

    x = np.arange(len(order))
    width = 0.36

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True,
                                   gridspec_kw={'height_ratios': [1, 1],
                                                'hspace': 0.18})

    for i, score in enumerate(['ePRS_Comp', 'PRS_G']):
        sub = d[d.score == score].set_index('label').loc[order]
        off = (i - 0.5) * width
        ax1.bar(x + off, sub[f'{tag}_ppv'] * 100, width * 0.92,
                color=SCORE_COLOURS[score], label=SCORE_LABELS[score],
                edgecolor='white', linewidth=0.6, zorder=3)
        ax2.bar(x + off, sub[f'{tag}_cases_caught'], width * 0.92,
                color=SCORE_COLOURS[score], edgecolor='white', linewidth=0.6,
                zorder=3)
        for xi, (p, c) in enumerate(zip(sub[f'{tag}_ppv'] * 100,
                                        sub[f'{tag}_cases_caught'])):
            ax1.text(xi + off, p, f'{p:.1f}', ha='center', va='bottom',
                     fontsize=7.5, color='#333333')
            ax2.text(xi + off, c, f'{int(c)}', ha='center', va='bottom',
                     fontsize=7.5, color='#333333')

    # Case rate in each subgroup: what flagging at random would return.
    base = d[d.score == 'ePRS_Comp'].set_index('label').loc[order, 'base_rate'] * 100
    for xi, b in enumerate(base):
        ax1.plot([xi - width, xi + width], [b, b], color='#444444', lw=1.4,
                 ls='--', zorder=4,
                 label='Case rate in subgroup' if xi == 0 else None)

    ax1.set_ylabel('Precision at the top 20 %\n(% of flagged who are cases)',
                   fontsize=9.5)
    ax1.set_title('A  Yield is concentrated above the subgroup case rate',
                  fontsize=10.5, fontweight='bold', loc='left', pad=6)
    ax1.legend(fontsize=8.5, frameon=True, edgecolor='#cccccc', loc='upper right')

    ax2.set_ylabel('Cases captured at the\ntop 20 % of each subgroup', fontsize=9.5)
    ax2.set_title('B  Number of cases that represents', fontsize=10.5,
                  fontweight='bold', loc='left', pad=6)

    n = d[d.score == 'ePRS_Comp'].set_index('label').loc[order]
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{lab}\nn={int(r.n):,}; {int(r.n_cases):,} cases'
                         for lab, r in zip(order, n.itertuples())], fontsize=8)

    for ax in (ax1, ax2):
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.yaxis.grid(True, color='#e8e8e8', zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(bottom=False)

    fig.text(0.01, -0.02,
             'Both scores flag the top 20 % of each subgroup, so they are compared at the same '
             'classification budget and neither benefits from flagging more people.\n'
             'Subgroups are ordered by case rate. "Not flagged" means the measure did not exceed '
             'its clinical threshold; "clinical model not flagged" uses age, sex and all three\n'
             'measures jointly rather than one cut-point at a time.',
             fontsize=7.5, color='#555555', va='top')

    fig.subplots_adjust(left=0.13, right=0.98, top=0.94, bottom=0.14)
    for ext in ('png', 'pdf'):
        out = os.path.join(args.out_dir, f'clinically_missed_subgroups.{ext}')
        fig.savefig(out, dpi=180, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)


if __name__ == '__main__':
    main()
