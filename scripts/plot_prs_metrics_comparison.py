#!/usr/bin/env python3
"""
plot_prs_metrics_comparison.py

Comprehensive PRS model comparison visualization.
Creates a 4-panel figure:
  - Panel A: Metrics heatmap across all pairwise comparisons
  - Panel B: Per-model performance (AUC, Recall, Precision)
  - Panel C: NRI diverging bar chart (who benefits from which model)
  - Panel D: Distinctiveness vs performance scatter (Jaccard vs AUC gap, sized by |NRI|)

Usage:
    from plot_prs_metrics_comparison import plot_prs_metrics_comparison

    plot_prs_metrics_comparison(
        mcnemar_path='McNemarStatsTestsAcrossPRSCalculations_refactored.csv',
        perf_path='model_recall_precision_improvement.csv',
        output_path='./figures'
    )

Or standalone:
    python plot_prs_metrics_comparison.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from typing import Optional
# Hatch legend entries
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib

import matplotlib.ticker as mticker

# ─── Pipeline color scheme ────────────────────────────────────────────────────
COHORT_COLORS = {
    'main':     '#E69F00',   # Orange
    'epi':      '#56B4E9',   # Sky blue
    'epi+main': '#CC79A7',   # Pinkish purple
    'cardio':   '#009E73',   # Bluish green
    'all':      '#F0E442',   # Yellow
    'combined': '#D55E00',   # Vermillion
    # Product variants (same as base)
    'main_product':     '#E69F00',
    'epi_product':      '#56B4E9',
    'epi+main_product': '#CC79A7',
    'cardio_product':   '#009E73',
    'all_product':      '#F0E442',
    # Summed variants (darker shades)
    'main_summed':      '#C77D00',  # Darker orange
    'epi_summed':       '#0073A8',  # Darker blue
    'epi+main_summed':  '#9F5580',  # Darker purple
    'cardio_summed':    '#006B52',  # Darker green
    'all_summed':       '#C7B800',  # Darker yellow
}

COHORT_MARKERS = {
    'main':     'o',
    'epi':      's',
    'epi+main': 'p',
    'cardio':   '^',
    'all':      'x'
}

# ─── Display labels for the performance/complementarity figure ────────────────
# Maps short model names to LaTeX-formatted PRS macro labels.
COHORT_DISPLAY_LABELS: dict = {
    'main':              r'$\mathrm{PRS}_G$',
    'epi':               r'$\mathrm{PRS}_{GxG}$',
    'epi+main':          r'$\mathrm{PRS}_{G+G}$',
    'cardio':            r'$\mathrm{PRS}_{GxGxE}$',
    'all':               r'$\mathrm{PRS}_{all}$',
    'combined':          r'$\mathrm{PRS}_{combined}$',
    # Explicit suffix variants
    'main_product':      r'$\mathrm{PRS}_G$',
    'epi_product':       r'$\mathrm{PRS}_{GxG}$',
    'epi+main_product':  r'$\mathrm{PRS}_{G+G}$',
    'cardio_product':    r'$\mathrm{PRS}_{GxGxE}$',
    'main_summed':       r'$\mathrm{PRS}_G$',
    'epi_summed':        r'$\mathrm{PRS}_{G+G}$',
    'epi+main_summed':   r'$\mathrm{PRS}_{G+G+G}$',
    'cardio_summed':     r'$\mathrm{PRS}_{G+GxE}$',
}

# Left-to-right display order for the performance figure (basic → complex)
PERF_MODEL_ORDER: list = ['main', 'epi+main', 'epi', 'cardio', 'all', 'combined']

# ─── Comparison display labels ────────────────────────────────────────────────
COMP_LABELS = {
    'main v epi':        'Main vs Epi',
    'main v cardio':     'Main vs Cardio',
    'main v epi+main':   'Main vs Epi+Main',
    'epi v cardio':      'Epi vs Cardio',
    'epi v epi+main':    'Epi vs Epi+Main',
    'cardio v epi+main': 'Cardio vs Epi+Main',
}

# Models in display order (base order for backward compatibility)
MODEL_ORDER = ['cardio', 'epi+main', 'epi', 'main']


def get_model_order(models_available: list) -> list:
    """
    Build display order with priority: base models first, then product, then summed.
    
    Order within each group: cardio, epi+main, epi, main
    """
    base_order = ['cardio', 'epi+main', 'epi', 'main']
    
    result = []
    # First add base models if present
    for base in base_order:
        if base in models_available:
            result.append(base)
    
    # Then add product variants
    for base in base_order:
        product = f'{base}_product'
        if product in models_available:
            result.append(product)
    
    # Then add summed variants
    for base in base_order:
        summed = f'{base}_summed'
        if summed in models_available:
            result.append(summed)
    
    # Add any remaining models not in the standard order
    for model in models_available:
        if model not in result:
            result.append(model)
    
    return result


# ─── Data extraction helpers ──────────────────────────────────────────────────

def _extract_stat(mcnemar: pd.DataFrame, comparison: str, test: str,
                  field: str = 'statistic') -> float:
    row = mcnemar[(mcnemar['comparison'] == comparison) & (mcnemar['test'] == test)]
    if row.empty:
        return np.nan
    return float(row[field].values[0])


def build_metrics_table(mcnemar: pd.DataFrame,
                        comparisons: list) -> pd.DataFrame:
    """
    Build a (comparisons x metrics) DataFrame for the heatmap.

    Metrics included:
        Cohen Kappa       – agreement (0=none, 1=perfect)
        Jaccard           – case-set overlap (0=distinct, 1=identical)
        Discordance       – proportion classified differently
        Pearson r         – PRS score correlation
        AUC diff (Δ)      – AUC(model2) − AUC(model1): + means model2 better
        NRI               – net reclassification improvement: + means model2 better
        NRI events        – NRI for true cases only
    """
    records = []
    for comp in comparisons:
        m1, m2 = [x.strip() for x in comp.split(' v ')]
        auc1 = _extract_stat(mcnemar, comp, 'DeLong', 'auc1')
        auc2 = _extract_stat(mcnemar, comp, 'DeLong', 'auc2')
        auc_diff = auc2 - auc1 if not np.isnan(auc1) else np.nan

        records.append({
            'comparison': comp,
            'label':      COMP_LABELS.get(comp, comp),
            'model1':     m1,
            'model2':     m2,
            'Kappa':              _extract_stat(mcnemar, comp, 'Cohen_Kappa'),
            'Jaccard':            _extract_stat(mcnemar, comp, 'Jaccard_Overlap'),
            'Discordance':        _extract_stat(mcnemar, comp, 'discordance_proportion'),
            'Pearson r':          _extract_stat(mcnemar, comp, 'pearson_r'),
            'AUC Δ (m2−m1)':      auc_diff,
            'NRI (overall)':      _extract_stat(mcnemar, comp, 'NRI'),
            'NRI (events)':       _extract_stat(mcnemar, comp, 'NRI', 'nri_events'),
        })
    return pd.DataFrame(records)


def filter_distinct_comparisons(
    metrics_df: pd.DataFrame,
    kappa_threshold: float = 0.4,
    priority_order: Optional[list] = None,
) -> pd.DataFrame:
    """
    Keep only comparisons between statistically distinct models.

    For any model with Kappa > kappa_threshold against multiple partners,
    retain only the single comparison with the highest-priority partner.

    Default priority: main > epi > cardio > epi+main > all > combined
    """
    if priority_order is None:
        priority_order = ['main', 'epi', 'cardio', 'epi+main', 'all', 'combined']
    priority_map = {m: i for i, m in enumerate(priority_order)}

    # Flag redundant pairs
    redundant_pairs = metrics_df[metrics_df['Kappa'] > kappa_threshold].copy()

    # Collect all models that appear in at least one redundant pair
    redundant_models = set(redundant_pairs['model1'].tolist() +
                           redundant_pairs['model2'].tolist())

    # For each redundant model, pick the single best-priority comparison to keep
    keep_comps: set = set()
    for model in redundant_models:
        involved = redundant_pairs[
            (redundant_pairs['model1'] == model) | (redundant_pairs['model2'] == model)
        ]
        # Pair each row with the other model's priority
        best_comp, best_pri = None, 999
        for _, row in involved.iterrows():
            partner = row['model2'] if row['model1'] == model else row['model1']
            pri = priority_map.get(partner, 999)
            if pri < best_pri:
                best_pri, best_comp = pri, row['comparison']
        if best_comp:
            keep_comps.add(best_comp)

    # Build final filtered list
    result = []
    for _, row in metrics_df.iterrows():
        if row['Kappa'] <= kappa_threshold:
            result.append(row)          # distinct pair — always keep
        elif row['comparison'] in keep_comps:
            result.append(row)          # redundant but chosen representative
        # else: skip silently

    return pd.DataFrame(result).reset_index(drop=True)


def build_model_perf(perf: pd.DataFrame) -> pd.DataFrame:
    """
    Clean up the performance dataframe, mapping to model short names.
    
    Handles both base models and suffixed variants (_product, _summed).
    Maps 'scaled_prs_cardio_product' -> 'cardio_product', etc.
    """
    df = perf.copy()
    
    # The cohort column contains strings like 'scaled_prs_cardio_product'
    # We want to extract everything after 'scaled_prs_'
    def extract_model(cohort_name):
        if pd.isna(cohort_name):
            return None
        if isinstance(cohort_name, str) and cohort_name.startswith('scaled_prs_'):
            return cohort_name.replace('scaled_prs_', '')
        return None
    
    df['model'] = df['cohort'].apply(extract_model)
    df = df.dropna(subset=['model'])
    
    # Keep one row per model (take first occurrence)
    df = df.groupby('model', sort=False).first()
    return df


# ─── Main plot function ───────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# Performance & Complementarity helpers
# ─────────────────────────────────────────────────────────────────────────────

def _hanley_mcneil_ci(auc: float, n_pos: int, n_neg: int,
                      z: float = 1.96) -> tuple[float, float]:
    """
    Compute a 95% CI for a single AUC using the Hanley–McNeil formula.

    Parameters
    ----------
    auc   : point estimate of AUC
    n_pos : number of positive cases in the validation set
    n_neg : number of negative controls in the validation set
    z     : z-score for desired CI width (default 1.96 → 95%)

    Returns
    -------
    (auc_lo, auc_hi)
    """
    if n_pos <= 0 or n_neg <= 0 or np.isnan(auc):
        return np.nan, np.nan
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc ** 2 / (1.0 + auc)
    se = np.sqrt(
        (auc * (1 - auc)
         + (n_pos - 1) * (q1 - auc ** 2)
         + (n_neg - 1) * (q2 - auc ** 2))
        / (n_pos * n_neg)
    )
    return max(0.0, auc - z * se), min(1.0, auc + z * se)


def build_performance_complementarity_table(
    mcnemar: pd.DataFrame,
    perf: pd.DataFrame,
    models_to_keep: Optional[list] = None,
) -> pd.DataFrame:
    """
    Build a per-model summary table for the performance & complementarity figure.

    Computes for every retained model:
      - AUC with Hanley–McNeil 95% CI
      - Cases identified in the top 20% risk bin
      - % of all cases identified (recall as fraction of total positives)
      - Exclusive cases (unique to this model)
      - Exclusive % of all cases
      - Precision and Recall (0–100 scale for plotting)
      - Display label (from COHORT_DISPLAY_LABELS)
      - Plot colour (from COHORT_COLORS)

    Parameters
    ----------
    mcnemar       : DataFrame from McNemarStatsTestsAcrossPRSCalculations_refactored.csv
    perf          : DataFrame from model_recall_precision_improvement.csv
    models_to_keep: list of short model names (e.g. ['main','epi','cardio']).
                    If None, all models present in perf are used.

    Returns
    -------
    pd.DataFrame with one row per model, columns:
        model, label, colour,
        auc, auc_lo, auc_hi,
        cases, total_cases, pct,
        exclusive, excl_pct,
        precision_pct, recall_pct
    """
    # ── 1. Strip 'scaled_prs_' prefix, map to short model name ───────────────
    def _to_short(cohort_name: str) -> Optional[str]:
        if pd.isna(cohort_name):
            return None
        s = str(cohort_name)
        return s.replace('scaled_prs_', '') if s.startswith('scaled_prs_') else s

    perf_clean = perf.copy()
    perf_clean['model'] = perf_clean['cohort'].apply(_to_short)
    perf_clean = perf_clean.dropna(subset=['model'])
    # One row per model — keep first occurrence (usually validation set)
    perf_clean = perf_clean.groupby('model', sort=False).first().reset_index()

    # ── 2. Optionally filter to kept models ───────────────────────────────────
    if models_to_keep is not None:
        perf_clean = perf_clean[perf_clean['model'].isin(models_to_keep)]

    # ── 3. Extract AUC from DeLong rows in mcnemar ────────────────────────────
    auc_map: dict[str, float] = {}
    delong = mcnemar[mcnemar['test'] == 'DeLong']
    for _, row in delong.iterrows():
        parts = row['comparison'].split(' v ')
        if len(parts) != 2:
            continue
        m1, m2 = parts[0].strip(), parts[1].strip()
        if not pd.isna(row.get('auc1')) and m1 not in auc_map:
            auc_map[m1] = float(row['auc1'])
        if not pd.isna(row.get('auc2')) and m2 not in auc_map:
            auc_map[m2] = float(row['auc2'])

    # ── 4. Build output records ───────────────────────────────────────────────
    records = []
    for _, pr in perf_clean.iterrows():
        model = pr['model']

        # Case counts (from threshold column at top-20% cut)
        tp   = int(pr.get('true_pos',  0))
        fn   = int(pr.get('false_neg', 0))
        fp   = int(pr.get('false_pos', 0))
        tn   = int(pr.get('true_neg',  0))
        n_pos = tp + fn          # total positives in validation set
        n_neg = fp + tn          # total negatives

        cases     = int(pr.get('cohort_high_risk_cases', tp))
        exclusive = int(pr.get('unique_cases', 0))
        pct       = (cases     / n_pos * 100) if n_pos > 0 else np.nan
        excl_pct  = (exclusive / n_pos * 100) if n_pos > 0 else np.nan

        # AUC + CI
        auc    = auc_map.get(model, np.nan)
        auc_lo, auc_hi = _hanley_mcneil_ci(auc, n_pos, n_neg)

        records.append({
            'model':        model,
            'label':        COHORT_DISPLAY_LABELS.get(model, model.upper()),
            'colour':       COHORT_COLORS.get(model, '#888888'),
            'auc':          auc,
            'auc_lo':       auc_lo,
            'auc_hi':       auc_hi,
            'cases':        cases,
            'total_cases':  n_pos,
            'pct':          pct,
            'exclusive':    exclusive,
            'excl_pct':     excl_pct,
            'precision_pct': float(pr.get('precision', np.nan)) * 100,
            'recall_pct':    float(pr.get('recall',    np.nan)) * 100,
        })

    result = pd.DataFrame(records)

    # ── 5. Sort by PERF_MODEL_ORDER (unknowns appended at end) ────────────────
    order_map = {m: i for i, m in enumerate(PERF_MODEL_ORDER)}
    result['_sort'] = result['model'].map(lambda m: order_map.get(m, 999))
    result = result.sort_values('_sort').drop(columns='_sort').reset_index(drop=True)

    return result


def plot_performance_complementarity(
    mcnemar_path: str,
    perf_path: str,
    output_path: str = '.',
    models_to_keep: Optional[list] = None,
    phenotype_label: str = '',
    multi_phenotype: Optional[list] = None,
    figsize: tuple = (18, 7.2),
    dpi: int = 300,
    verbose: bool = True,
) -> plt.Figure:
    """
    Two-panel performance and complementarity figure.

    Panel A  — Bars: % of all cases identified in top-20% risk bin, with absolute
               case counts annotated on each bar.  Exclusive cases shown as a
               hatched overlay.  AUC overlaid via a twin right y-axis (diamond
               markers with Hanley–McNeil 95% CI whiskers).

    Panel B  — Grouped bar chart: precision (lighter fill) and recall (hatched)
               side-by-side per model, values annotated above each bar.

    Parameters
    ----------
    mcnemar_path     : Path to McNemarStatsTestsAcrossPRSCalculations_refactored.csv.
                       Ignored when multi_phenotype is provided.
    perf_path        : Path to model_recall_precision_improvement.csv.
                       Ignored when multi_phenotype is provided.
    output_path      : Directory to save the figure.
    models_to_keep   : Short model names to include (e.g. ['main','epi','cardio']).
                       If None, all models in perf are included.
    phenotype_label  : Disease / phenotype name shown below the bars (single-
                       phenotype mode). Defaults to empty string.
    multi_phenotype  : Optional list of dicts, one per phenotype, each with keys:
                           label          : str  — phenotype display name
                           mcnemar_path   : str  — path to McNemar CSV
                           perf_path      : str  — path to perf CSV
                           models_to_keep : list — models to include (or None)
                       When provided, mcnemar_path / perf_path / models_to_keep /
                       phenotype_label at the top level are ignored.
    figsize          : Figure dimensions (width, height).
    dpi              : Output resolution.
    verbose          : Print save path.

    Returns
    -------
    plt.Figure
    """
    # ── Assemble per-phenotype data ───────────────────────────────────────────
    if multi_phenotype is not None:
        groups: list[tuple[str, pd.DataFrame]] = []
        for entry in multi_phenotype:
            mc  = pd.read_csv(entry['mcnemar_path'])
            pf  = pd.read_csv(entry['perf_path'])
            mtk = entry.get('models_to_keep', None)
            tbl = build_performance_complementarity_table(mc, pf, mtk)
            groups.append((entry['label'], tbl))
    else:
        mc  = pd.read_csv(mcnemar_path)
        pf  = pd.read_csv(perf_path)
        tbl = build_performance_complementarity_table(mc, pf, models_to_keep)
        groups = [(phenotype_label, tbl)]

    # ── Figure + axes ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=figsize, facecolor='white')
    import matplotlib.gridspec as _gs
    grid = _gs.GridSpec(
        1, 2, figure=fig,
        left=0.06, right=0.97,
        bottom=0.18, top=0.90,
        wspace=0.38,
        width_ratios=[1.45, 1.0],
    )
    ax_bar = fig.add_subplot(grid[0])
    ax_auc = ax_bar.twinx()        # right axis — AUC
    ax_pr  = fig.add_subplot(grid[1])

    for ax in [ax_bar, ax_pr]:
        ax.set_facecolor('white')
        ax.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_pr.spines['right'].set_visible(False)
    ax_pr.spines['top'].set_visible(False)

    BAR_W       = 0.72
    GRP_W       = 0.38    # width of each P/R bar within a pair
    GRP_GAP     = 0.04    # gap between precision and recall bars
    MODEL_GAP   = 0.55    # gap between model groups in Panel B
    DISEASE_GAP = 1.3     # extra gap between disease groups in Panel B
    DISEASE_SEP = 1.0     # extra x-gap between disease groups in Panel A

    x_a   = 0   # running x position for Panel A
    x_b   = 0   # running x position for Panel B
    a_xticks, a_xlabels, a_disease_xmid = [], [], {}
    b_xticks, b_xlabels, b_disease_xmid = [], [], {}

    for pheno_label, tbl in groups:
        xs_a, xs_b = [], []

        for _, d in tbl.iterrows():
            col = d['colour']
            alp = 1.0   # single-phenotype always full opacity; caller controls via colour

            # ── Panel A bars ─────────────────────────────────────────────────
            # Total % (light fill)
            ax_bar.bar(x_a, d['pct'], width=BAR_W,
                       color=col, alpha=0.38,
                       edgecolor=col, linewidth=0.9, zorder=2)
            # Exclusive % (hatched overlay)
            ax_bar.bar(x_a, d['excl_pct'], width=BAR_W,
                       color=col, alpha=alp,
                       hatch='////', edgecolor='white', linewidth=0.5, zorder=3)
            # Absolute n annotated above bar
            ax_bar.text(x_a, d['pct'] + 0.6,
                        f'n={int(d["cases"]):,}',
                        ha='center', va='bottom', fontsize=5.8,
                        color=col, fontweight='bold')

            # ── Panel A AUC overlay (twin axis) ──────────────────────────────
            # CI whisker
            ax_auc.plot([x_a, x_a], [d['auc_lo'], d['auc_hi']],
              color=col, lw=1.4, alpha=0.65, zorder=5)
            ax_auc.plot([x_a - 0.18, x_a + 0.18], [d['auc_lo']] * 2,
              color=col, lw=1.0, alpha=0.55, zorder=5)
            ax_auc.plot([x_a - 0.18, x_a + 0.18], [d['auc_hi']] * 2,
              color=col, lw=1.0, alpha=0.55, zorder=5)
            # Diamond marker
            ax_auc.scatter(x_a, d['auc'],
              color=col, s=70, marker='D',
              alpha=alp, edgecolors='white', linewidth=1.0, zorder=6)
            # Value label right of marker
            ax_auc.text(x_a, d['auc'] - .03,
              f'{d["auc"]:.2f}',
              va='center', ha='left', fontsize=6.2,
              color=col, style='italic')
  
            a_xticks.append(x_a)
            a_xlabels.append(d['label'])
            xs_a.append(x_a)
            x_a += 1

            # ── Panel B precision / recall bars ──────────────────────────────
            xp = x_b
            xr = x_b + GRP_W + GRP_GAP
            ax_pr.bar(xp, d['precision_pct'], width=GRP_W,
                      color=col, alpha=0.50,
                      edgecolor=col, linewidth=0.8, zorder=2)
            ax_pr.bar(xr, d['recall_pct'], width=GRP_W,
                      color=col, alpha=alp,
                      hatch='XXX', edgecolor='white', linewidth=0.5, zorder=3)
            ax_pr.text(xp, d['precision_pct'] + 0.7,
                       f'{d["precision_pct"]:.0f}%',
                       ha='center', va='bottom', fontsize=6.5,
                       color=col, fontweight='bold')
            ax_pr.text(xr, d['recall_pct'] + 0.7,
                       f'{d["recall_pct"]:.0f}%',
                       ha='center', va='bottom', fontsize=6.5,
                       color=col, fontweight='bold')

            tick_b = x_b + GRP_W + GRP_GAP / 2
            b_xticks.append(tick_b)
            b_xlabels.append(d['label'])
            xs_b.append(tick_b)
            x_b += GRP_W * 2 + GRP_GAP + MODEL_GAP

        a_disease_xmid[pheno_label] = np.mean(xs_a)
        b_disease_xmid[pheno_label] = np.mean(xs_b)
        x_a += DISEASE_SEP
        x_b += DISEASE_GAP

    # ── Panel A axis formatting ───────────────────────────────────────────────
    ax_bar.set_ylim(0, 82)
    ax_bar.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))
    ax_bar.set_ylabel(
        'Cases identified in top 20% risk\n(% of all cases in validation set)',
        fontsize=9.0)
    ax_bar.set_xticks(a_xticks)
    ax_bar.set_xticklabels(a_xlabels, fontsize=8.0, rotation=35, ha='right')
    ax_bar.set_title('A', loc='left', fontsize=12, fontweight='bold')

    ax_auc.set_ylim(0.45, 0.98)
    ax_auc.set_ylabel('AUC (95% CI)', fontsize=9.0, color='#444444')
    ax_auc.tick_params(axis='y', labelcolor='#444444')
    ax_auc.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{v:.2f}'))
    ax_auc.spines['top'].set_visible(False)
    ax_auc.axhline(0.5, color='#cccccc', lw=0.7, ls=':', zorder=0)

    # Disease group labels below Panel A ticks
    if len(groups) > 1:
        for lbl, xmid in a_disease_xmid.items():
            ax_bar.annotate(
                lbl,
                xy=(xmid, 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -52), textcoords='offset points',
                ha='center', va='top',
                fontsize=9.0, fontweight='bold', color='#333333',
                annotation_clip=False)
        # Vertical separator between disease groups
        sep_x = len(groups[0][1]) - 0.5 + DISEASE_SEP / 2
        ax_bar.axvline(sep_x, color='#cccccc', lw=1.0, ls='--', zorder=0)
    elif groups[0][0]:   # single phenotype with a label
        ax_bar.set_xlabel(groups[0][0], fontsize=9.5, fontweight='bold',
                          color='#333333')

    # Legends for Panel A
    ax_bar.legend(
        handles=[
            mpatches.Patch(facecolor='#aaaaaa', alpha=0.38, edgecolor='#aaaaaa',
                           label='Total cases identified'),
            mpatches.Patch(facecolor='#aaaaaa', alpha=0.90, hatch='////',
                           edgecolor='white', label='Exclusive to model'),
        ],
        fontsize=7.5, loc='upper left', framealpha=0.88, edgecolor='#cccccc')

    ax_auc.legend(
        handles=[
            plt.Line2D([0], [0], marker='D', color='#555555', lw=0,
                       ms=7, label='AUC (95% CI)'),
        ],
        fontsize=7.5, loc='upper right', framealpha=0.88, edgecolor='#cccccc')

    # ── Panel B axis formatting ───────────────────────────────────────────────
    ax_pr.set_xticks(b_xticks)
    ax_pr.set_xticklabels(b_xlabels, fontsize=8.0, rotation=35, ha='right')
    ax_pr.set_ylabel('Precision / Recall (%)', fontsize=9.0)
    ax_pr.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{v:.0f}%'))
    ax_pr.set_ylim(0, 97)
    ax_pr.set_title('B', loc='left', fontsize=12, fontweight='bold')

    if len(groups) > 1:
        for lbl, xmid in b_disease_xmid.items():
            ax_pr.annotate(
                lbl,
                xy=(xmid, 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -52), textcoords='offset points',
                ha='center', va='top',
                fontsize=9.0, fontweight='bold', color='#333333',
                annotation_clip=False)
        n_first = len(groups[0][1])
        sep_xb = b_xticks[n_first - 1] + GRP_W + DISEASE_GAP * 0.5
        ax_pr.axvline(sep_xb, color='#cccccc', lw=1.0, ls='--', zorder=0)
    elif groups[0][0]:
        ax_pr.set_xlabel(groups[0][0], fontsize=9.5, fontweight='bold',
                         color='#333333')

    ax_pr.legend(
        handles=[
            mpatches.Patch(facecolor='#aaaaaa', alpha=0.50, edgecolor='#aaaaaa',
                           label='Precision'),
            mpatches.Patch(facecolor='#aaaaaa', alpha=0.95, hatch='XXX',
                           edgecolor='white', label='Recall'),
        ],
        fontsize=7.5, loc='upper left', framealpha=0.88, edgecolor='#cccccc')

    # ── Shared model colour legend (bottom of figure) ─────────────────────────
    all_models = []
    seen = set()
    for _, tbl in groups:
        for _, d in tbl.iterrows():
            if d['model'] not in seen:
                all_models.append(d)
                seen.add(d['model'])

    fig.legend(
        handles=[
            mpatches.Patch(color=d['colour'], label=d['label'])
            for d in all_models
        ],
        loc='lower center', ncol=len(all_models),
        fontsize=8.2, framealpha=0.92, edgecolor='#cccccc',
        bbox_to_anchor=(0.5, 0.005))

    # ── Figure title ─────────────────────────────────────────────────────────
    title = 'ePRS Model Performance and Complementarity'
    if len(groups) > 1:
        title += ' — ' + ' and '.join(lbl for lbl, _ in groups if lbl)
    elif groups[0][0]:
        title += f' — {groups[0][0]}'
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.975)

    # ── Save ─────────────────────────────────────────────────────────────────
    os.makedirs(output_path, exist_ok=True)
    out_file = os.path.join(output_path, 'prs_performance_complementarity.png')
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    out_pdf = os.path.join(output_path, 'prs_performance_complementarity.pdf')
    fig.savefig(out_pdf, dpi=dpi, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    if verbose:
        print(f'  ✓ Saved: {out_file}')

    return fig


# ─── Main plot function ───────────────────────────────────────────────────────

def plot_prs_metrics_comparison(
    mcnemar_path: str,
    perf_path: str,
    output_path: str = '.',
    models_to_keep: Optional[list] = None,
    auc_map: Optional[dict] = None,
    exclude_comparisons: Optional[list] = None,
    figsize: tuple = (20, 15),
    dpi: int = 300,
    verbose: bool = True,
) -> plt.Figure:
    """
    Create a 3-panel PRS model comparison metrics figure.

    Parameters
    ----------
    mcnemar_path : str
        Path to McNemarStatsTests CSV file.
    perf_path : str
        Path to model_recall_precision_improvement CSV file.
    output_path : str
        Directory to save the figure.
    models_to_keep : list, optional
        Models confirmed statistically distinct (from quick_filter_models /
        filter_redundant_models). When provided, Panels B and C are filtered
        to these models only. The heatmap still shows one representative row
        per redundant model so the redundancy decision is visible.
        If None, redundant models are detected internally via Kappa threshold.
    auc_map : dict, optional
        Override AUC values per model e.g. {'cardio': 0.825, 'main': 0.591}.
        If None, AUCs are extracted from the DeLong rows.
    exclude_comparisons : list, optional
        Comparisons to exclude (default: excludes 'all v combined').
    figsize : tuple
        Figure dimensions (width, height).
    dpi : int
        Output resolution.
    verbose : bool
        Print progress.

    Returns
    -------
    plt.Figure
    """
    # ── Load data ─────────────────────────────────────────────────────────────
    mcnemar = pd.read_csv(mcnemar_path)
    perf    = pd.read_csv(perf_path)

    if exclude_comparisons is None:
        exclude_comparisons = ['all v combined']

    all_comparisons = [c for c in mcnemar['comparison'].unique()
                       if c not in exclude_comparisons]

    # ── Filter to distinct comparisons only ───────────────────────────────────
    # A model is "redundant" if it has Kappa >= threshold with ANY partner.
    # For such models, ALL their comparisons are dropped except the ONE
    # with the highest-priority partner (main > epi > cardio).
    PRIORITY     = {'main': 0, 'epi': 1, 'cardio': 2, 'epi+main': 3,
                    'all': 4, 'combined': 5}
    KAPPA_THRESH = 0.4

    # Step 1: identify redundant models (high Kappa with at least one partner)
    redundant_models = set()
    for comp in all_comparisons:
        kappa = _extract_stat(mcnemar, comp, 'Cohen_Kappa')
        if not np.isnan(kappa) and kappa >= KAPPA_THRESH:
            m1, m2 = [x.strip() for x in comp.split(' v ')]
            # The lower-priority model is the redundant one
            red = m2 if PRIORITY.get(m2, 99) > PRIORITY.get(m1, 99) else m1
            redundant_models.add(red)

    # If caller already knows which models to keep, derive redundant_models from that
    # instead of relying solely on internal Kappa detection.
    if models_to_keep is not None:
        all_model_names = set()
        for comp in all_comparisons:
            for m in comp.split(' v '):
                all_model_names.add(m.strip())
        redundant_models = all_model_names - set(models_to_keep)

    if verbose:
        print(f"  Redundant models identified: {redundant_models}")
        if models_to_keep is not None:
            print(f"  (overridden by models_to_keep={models_to_keep})")

    # Step 2: for each redundant model choose best-priority representative
    kept_redundant = {}   # redundant_model -> comparison to keep
    for comp in all_comparisons:
        m1, m2 = [x.strip() for x in comp.split(' v ')]
        for red in redundant_models:
            if red in (m1, m2):
                partner = m2 if m1 == red else m1
                pri = PRIORITY.get(partner, 99)
                if red not in kept_redundant or pri < PRIORITY.get(
                        (kept_redundant[red].split(' v ')[1].strip()
                         if kept_redundant[red].split(' v ')[0].strip() == red
                         else kept_redundant[red].split(' v ')[0].strip()), 99):
                    kept_redundant[red] = comp

    # Step 3: build final comparison list
    comparisons = []
    for comp in all_comparisons:
        m1, m2 = [x.strip() for x in comp.split(' v ')]
        involves_redundant = any(r in (m1, m2) for r in redundant_models)
        if not involves_redundant:
            comparisons.append(comp)
        elif comp in kept_redundant.values():
            comparisons.append(comp)

    if verbose:
        print(f"  Kept comparisons: {comparisons}")

    metrics_df  = build_metrics_table(mcnemar, comparisons)
    model_perf  = build_model_perf(perf)

    # AUC per model: derive from DeLong rows (first occurrence as model1 or model2)
    if auc_map is None:
        auc_map = {}
        delong = mcnemar[mcnemar['test'] == 'DeLong']
        for _, row in delong.iterrows():
            m1, m2 = [x.strip() for x in row['comparison'].split(' v ')]
            if not pd.isna(row['auc1']) and m1 not in auc_map:
                auc_map[m1] = float(row['auc1'])
            if not pd.isna(row['auc2']) and m2 not in auc_map:
                auc_map[m2] = float(row['auc2'])

    # ── Layout ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=figsize, facecolor='#F7F9FC')
    gs  = gridspec.GridSpec(
        2, 2,
        figure=fig,
        hspace=0.52,
        wspace=0.38,
        left=0.17, right=0.97,
        top=0.88, bottom=0.07
    )
    ax_heat   = fig.add_subplot(gs[0, :])   # Full top row: heatmap
    ax_perf   = fig.add_subplot(gs[1, 0])   # Bottom-left: model performance
    ax_nri    = fig.add_subplot(gs[1, 1])   # Bottom-right: NRI diverging bars

    axes = [ax_heat, ax_perf, ax_nri]
    for ax in axes:
        ax.set_facecolor('#FFFFFF')
        for spine in ax.spines.values():
            spine.set_edgecolor('#D0D0D0')
            spine.set_linewidth(0.8)

    # ─────────────────────────────────────────────────────────────────────────
    # PANEL A: Metrics Heatmap
    # ─────────────────────────────────────────────────────────────────────────
    metric_cols = ['Kappa', 'Jaccard', 'Discordance', 'Pearson r',
                   'AUC Δ (m2−m1)']
                
#   metric_cols = ['Kappa', 'Jaccard', 'Discordance', 'Pearson r',
#                  'AUC Δ (m2−m1)', 'NRI (overall)', 'NRI (events)']

    heatmap_data = metrics_df.set_index('label')[metric_cols]

    # Each metric has its own colour interpretation:
    # Kappa & Jaccard    – 0=good (distinct), 1=bad (redundant) → Oranges
    # Discordance        – higher = more distinct → Blues
    # Pearson r          – diverging around 0
    # AUC Δ & NRI cols   – diverging around 0; positive = model2 better
    colormaps = {
        'Kappa':            ('Oranges', 0.0,  1.0,  False),
        'Jaccard':          ('Oranges', 0.0,  1.0,  False),
        'Discordance':      ('Blues',   0.0,  1.0,  False),
        'Pearson r':        ('RdBu',   -0.6,  0.6,  True),
        'AUC Δ (m2−m1)':   ('RdBu',   -0.35, 0.35, True)
#       'NRI (overall)':   ('RdBu',   -0.5,  0.5,  True),
#       'NRI (events)':    ('RdBu',   -0.45, 0.45, True),
    }

    n_rows, n_cols = heatmap_data.shape
    col_w = 1.0 / n_cols
    row_h = 1.0 / n_rows

    # Draw each cell individually with per-column colouring
    for ci, metric in enumerate(metric_cols):
        cmap_name, vmin, vmax, diverging = colormaps[metric]
        cmap = plt.get_cmap(cmap_name)
        vals = heatmap_data[metric].values

        for ri, val in enumerate(vals):
            if np.isnan(val):
                color = '#EEEEEE'
                text_color = '#AAAAAA'
                display = 'N/A'
            else:
                if diverging:
                    norm_val = (val - vmin) / (vmax - vmin)
                else:
                    norm_val = (val - vmin) / (vmax - vmin)
                norm_val = np.clip(norm_val, 0, 1)
                color = cmap(norm_val)
                # Text contrast
                brightness = 0.299*color[0] + 0.587*color[1] + 0.114*color[2]
                text_color = 'white' if brightness < 0.55 else '#1A1A1A'
                display = f'{val:+.2f}' if diverging else f'{val:.2f}'

            x = ci / n_cols
            y = (n_rows - 1 - ri) / n_rows
            rect = mpatches.FancyBboxPatch(
                (x + 0.003, y + 0.015),
                col_w - 0.006, row_h - 0.025,
                boxstyle="round,pad=0.01",
                linewidth=0,
                facecolor=color,
                transform=ax_heat.transAxes,
                zorder=2
            )
            ax_heat.add_patch(rect)
            ax_heat.text(
                x + col_w / 2, y + row_h / 2,
                display,
                ha='center', va='center',
                fontsize=9.5, fontweight='bold',
                color=text_color,
                transform=ax_heat.transAxes,
                zorder=3
            )

    # Y-axis labels — render "M1 COLOR  vs  M2 COLOR" as separate text spans
    # rendered RIGHT → LEFT from x = -0.002 in axes coordinates
    CHAR_W = 0.005   # approx axes-coord width per uppercase character at font 10/Arial
    VS_W   = 0.022    # width of " vs " at font size 9
    GAP    = 0.002    # gap between text pieces

    for ri, label in enumerate(heatmap_data.index):
        y = (n_rows - 1 - ri) / n_rows + row_h / 2
        m1 = metrics_df.set_index('label').loc[label, 'model1']
        m2 = metrics_df.set_index('label').loc[label, 'model2']
        c1 = COHORT_COLORS.get(m1, '#1A1A1A')
        c2 = COHORT_COLORS.get(m2, '#1A1A1A')

        m1u, m2u = m1.upper(), m2.upper()

        # 1) M2 — right-aligned at x = -0.002
        ax_heat.text(-0.002, y, m2u,
                     ha='right', va='center', fontsize=8, fontweight='bold',
                     color=c2, transform=ax_heat.transAxes)

        # 2) " vs " — right edge lands just before M2 left edge
        x_vs = -0.002 - len(m2u) * CHAR_W - GAP
        ax_heat.text(x_vs, y, 'vs',
                     ha='right', va='center', fontsize=7, color='#777777',
                     transform=ax_heat.transAxes)

        # 3) M1 — right edge lands just before " vs " left edge
        x_m1 = x_vs - VS_W - GAP
        ax_heat.text(x_m1, y, m1u,
                     ha='right', va='center', fontsize=8, fontweight='bold',
                     color=c1, transform=ax_heat.transAxes)

    # X-axis labels (metrics)
    metric_display = {
        'Kappa':           "Cohen's\nKappa",
        'Jaccard':         'Jaccard\nOverlap',
        'Discordance':     'Discordance\nProportion',
        'Pearson r':       'Pearson r\n(PRS scores)',
        'AUC Δ (m2−m1)':  'AUC Δ\n(m2 − m1)'
#       'NRI (overall)':   'NRI\n(Overall)',
#       'NRI (events)':    'NRI\n(Events)',
    }
    for ci, metric in enumerate(metric_cols):
        x = ci / n_cols + col_w / 2
        ax_heat.text(
            x, 1.02,
            metric_display.get(metric, metric),
            ha='center', va='bottom',
            fontsize=9.5, fontweight='bold', color='#333333',
            transform=ax_heat.transAxes
        )

    ax_heat.set_xlim(0, 1)
    ax_heat.set_ylim(0, 1)
    ax_heat.axis('off')

    # Subtitle annotations for the heatmap legend
    legend_items = [
        mpatches.Patch(color=plt.get_cmap('Oranges')(0.7), label='Kappa / Jaccard  (darker = more overlap / redundant)'),
        mpatches.Patch(color=plt.get_cmap('Blues')(0.7),   label='Discordance  (darker = more distinct)'),
        mpatches.Patch(color=plt.get_cmap('RdBu')(0.15),   label='Pearson r / AUC \u0394 / NRI  (red = model 2 better)'),
        mpatches.Patch(color=plt.get_cmap('RdBu')(0.85),   label='(blue = model 1 better)'),
    ]
    ax_heat.legend(
        handles=legend_items,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.14),
        ncol=4,
        fontsize=8,
        frameon=True,
        framealpha=0.9,
        edgecolor='#CCCCCC',
        facecolor='#F7F9FC',
    )

    # No title for the heatmap panel

    # ─────────────────────────────────────────────────────────────────────────
    # PANEL B: Model Performance — distinct models only
    # ─────────────────────────────────────────────────────────────────────────
    # Get all available models and order them properly
    available_models = [m for m in model_perf.index if m not in redundant_models]
    models_ordered = get_model_order(available_models)
    colors_ordered = [COHORT_COLORS.get(m, '#888888') for m in models_ordered]

    # Metrics to show
    def _scalar(df, model, col):
        """Safely extract a single float from a potentially multi-row index."""
        if model not in df.index:
            return np.nan
        val = df.loc[model, col]
        return float(val.iloc[0]) if hasattr(val, 'iloc') else float(val)

    perf_metrics = {
        'AUC':       [auc_map.get(m, np.nan) for m in models_ordered],
        'Recall':    [_scalar(model_perf, m, 'recall')    for m in models_ordered],
        'Precision': [_scalar(model_perf, m, 'precision') for m in models_ordered],
    }

    x      = np.arange(len(models_ordered))
    n_bars = len(perf_metrics)
    w      = 0.22
    offsets = np.linspace(-(n_bars - 1) * w / 2, (n_bars - 1) * w / 2, n_bars)

    metric_styles = {
        'AUC':       dict(alpha=0.92, edgecolor='none',   hatch=''),
        'Recall':    dict(alpha=0.70, edgecolor='white',  hatch='//'),
        'Precision': dict(alpha=0.55, edgecolor='white',  hatch='..'),
    }

    for j, (metric_name, vals) in enumerate(perf_metrics.items()):
        style = metric_styles[metric_name]
        bars = ax_perf.bar(
            x + offsets[j], vals,
            width=w,
            color=colors_ordered,
            label=metric_name,
            **style,
            linewidth=0.5,
            zorder=3
        )
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax_perf.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.012,
                    f'{val:.2f}',
                    ha='center', va='bottom',
                    fontsize=7.5, color='#333333', fontweight='semibold'
                )


    legend_handles = [
        Patch(facecolor='#888888', alpha=0.92, label='AUC'),
        Patch(facecolor='#888888', alpha=0.70, hatch='//', edgecolor='white', label='Recall'),
        Patch(facecolor='#888888', alpha=0.55, hatch='..', edgecolor='white', label='Precision'),
    ]
    # Model colour swatches
    model_handles = [
        Patch(facecolor=COHORT_COLORS.get(m, '#888'), label=m.upper())
        for m in models_ordered
    ]
    ax_perf.legend(
        handles=legend_handles + model_handles,
        fontsize=7.5, ncol=2, framealpha=0.85,
        edgecolor='#DDDDDD', loc='upper right'
    )

    # Reference lines
    ax_perf.axhline(0.5,  color='#AAAAAA', lw=0.8, ls='--', zorder=1)
    ax_perf.axhline(0.75, color='#CCCCCC', lw=0.6, ls=':', zorder=1)

    ax_perf.set_xticks(x)
    ax_perf.set_xticklabels([m.upper() for m in models_ordered], fontsize=11, fontweight='bold')
    for tick, m in zip(ax_perf.get_xticklabels(), models_ordered):
        tick.set_color(COHORT_COLORS.get(m, '#333333'))
    ax_perf.set_ylabel('Metric Value', fontsize=8)
    ax_perf.set_ylim(0, 1.0)
    ax_perf.set_title(
        'Panel B  —  Per-Model Performance',
        fontsize=12, fontweight='bold', color='#1A1A1A', loc='left', pad=8
    )
    ax_perf.yaxis.grid(True, color='#EEEEEE', linewidth=0.7, zorder=0)
    ax_perf.set_axisbelow(True)

    # Annotate random-level threshold
    ax_perf.text(
        0.01, 0.51, '',
        transform=ax_perf.get_yaxis_transform(),
        fontsize=7, color='#AAAAAA', va='bottom'
    )

    # ─────────────────────────────────────────────────────────────────────────
    # PANEL C: NRI — distinct-model comparisons only
    # ─────────────────────────────────────────────────────────────────────────
    nri_df   = metrics_df[~metrics_df['model1'].isin(redundant_models) &
                          ~metrics_df['model2'].isin(redundant_models)].copy()
    nri_data = nri_df.set_index('label')[['NRI (overall)', 'NRI (events)',
                                           'model1', 'model2']].copy()
    nri_data = nri_data.sort_values('NRI (overall)')

    y_pos = np.arange(len(nri_data))
    bar_h = 0.35

    def nri_color(val):
        if np.isnan(val):
            return '#CCCCCC'
        return '#2166AC' if val > 0 else '#D6604D'

    # ── Pre-compute axis limits BEFORE drawing so labels stay inside ──────────
    # Padding is sized to fit the label text, not a % of the bar length.
    # "+0.00" is ~5 chars; at fontsize 8 in a ~10-inch panel that is roughly
    # 5% of the data span. We add a small gap between bar tip and text (1.5%).
    all_vals  = [v for v in list(nri_data['NRI (overall)'].dropna()) +
                              list(nri_data['NRI (events)'].dropna())]
    v_lo_data = min(min(all_vals, default=-0.1), 0)   # always include zero
    v_hi_data = max(max(all_vals, default=0.1),  0)

    # Rough axis span before padding — used to size gap and text allowance
    raw_span  = max(abs(v_lo_data) + abs(v_hi_data), 0.1)
    lbl_gap   = raw_span * 0.015          # gap between bar tip and label start
    lbl_width = raw_span * 0.055          # approximate width of "+0.00" text
    side_pad  = lbl_gap + lbl_width       # total space needed outside bar tip

    # Each side only gets padding if it actually has values; otherwise minimal
    pad_lo = side_pad if v_lo_data < 0 else raw_span * 0.03
    pad_hi = side_pad if v_hi_data > 0 else raw_span * 0.03
    x_lo   = v_lo_data - pad_lo
    x_hi   = v_hi_data + pad_hi
    ax_nri.set_xlim(x_lo, x_hi)

    # ── Draw bars + inline labels ─────────────────────────────────────────────
    for i, (label, row) in enumerate(nri_data.iterrows()):
        v_overall = row['NRI (overall)']
        v_events  = row['NRI (events)']
        if not np.isnan(v_overall):
            ax_nri.barh(y_pos[i] + bar_h / 2, v_overall, height=bar_h,
                        color=nri_color(v_overall), alpha=0.85,
                        label='NRI Overall' if i == 0 else '', zorder=3)
            ax_nri.text(
                v_overall + (lbl_gap if v_overall >= 0 else -lbl_gap),
                y_pos[i] + bar_h / 2,
                f'{v_overall:+.2f}',
                va='center', ha='left' if v_overall >= 0 else 'right',
                fontsize=8, fontweight='bold',
                color=nri_color(v_overall),
                clip_on=False,
            )
        if not np.isnan(v_events):
            ax_nri.barh(y_pos[i] - bar_h / 2, v_events, height=bar_h,
                        color=nri_color(v_events), alpha=0.50, hatch='////',
                        edgecolor='white', linewidth=0.3,
                        label='NRI Events' if i == 0 else '', zorder=3)

    ax_nri.axvline(0, color='#555555', lw=1.2, zorder=4)

    # ── Y-axis: two-tone M1 COLOR  vs  M2 COLOR labels ───────────────────────
    ax_nri.set_yticks(y_pos)
    ax_nri.set_yticklabels([''] * len(y_pos))

    for i, (label, row) in enumerate(nri_data.iterrows()):
        m1 = row['model1']
        m2 = row['model2']
        c1 = COHORT_COLORS.get(m1, '#333333')
        c2 = COHORT_COLORS.get(m2, '#333333')
        # Render right-to-left: M2 | vs | M1, right-aligned to left axis edge
        ax_nri.text(-0.02, y_pos[i], m2.upper(),
                    ha='right', va='center', fontsize=8, fontweight='bold',
                    color=c2, transform=ax_nri.get_yaxis_transform(), clip_on=False)
        ax_nri.text(-0.11, y_pos[i], 'vs',
                    ha='right', va='center', fontsize=7, color='#888888',
                    transform=ax_nri.get_yaxis_transform(), clip_on=False)
        ax_nri.text(-0.15, y_pos[i], m1.upper(),
                    ha='right', va='center', fontsize=8, fontweight='bold',
                    color=c1, transform=ax_nri.get_yaxis_transform(), clip_on=False)

    ax_nri.set_xlabel('NRI  (positive = model 2 improves on model 1)', fontsize=9)
    ax_nri.xaxis.grid(True, color='#EEEEEE', linewidth=0.7, zorder=0)
    ax_nri.set_axisbelow(True)
    ax_nri.set_title(
        'Panel C  —  Net Reclassification Improvement',
        fontsize=12, fontweight='bold', color='#1A1A1A', loc='left', pad=8
    )

    nri_legend = [
        mpatches.Patch(color='#2166AC', alpha=0.85, label='Positive NRI  (m2 better)'),
        mpatches.Patch(color='#D6604D', alpha=0.85, label='Negative NRI  (m1 better)'),
        mpatches.Patch(facecolor='#888888', alpha=0.50, hatch='////',
                       edgecolor='white', label='NRI events (hatch)'),
    ]
    ax_nri.legend(handles=nri_legend, fontsize=8, framealpha=0.85,
                  edgecolor='#DDDDDD', loc='lower right')

    # ─────────────────────────────────────────────────────────────────────────
    # Figure-level title and annotation
    # ─────────────────────────────────────────────────────────────────────────
    fig.suptitle(
        'PRS Model Comparison — Distinctiveness & Performance Metrics',
        fontsize=16, fontweight='bold', color='#1A1A1A', y=0.985
    )
    fig.text(
        0.5, 0.962,
        'Validation cohort  |  Threshold: top decile (bin > 8)  |  N = 1,569 main cases',
        ha='center', fontsize=10, color='#666666', style='italic'
    )

    # ─────────────────────────────────────────────────────────────────────────
    # Save
    # ─────────────────────────────────────────────────────────────────────────
    os.makedirs(output_path, exist_ok=True)
    out_file = os.path.join(output_path, 'prs_metrics_comparison.png')
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight', facecolor=fig.get_facecolor())

    if verbose:
        print(f'  ✓ Saved: {out_file}')

    return fig


# ─── CLI / standalone ─────────────────────────────────────────────────────────
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Create PRS metrics comparison figure'
    )
    parser.add_argument('--mcnemar', default='McNemarStatsTestsAcrossPRSCalculations_refactored.csv')
    parser.add_argument('--perf',    default='model_recall_precision_improvement.csv')
    parser.add_argument('--output',  default='./figures')
    parser.add_argument('--dpi',     type=int, default=300)
    args = parser.parse_args()

#   fig = plot_prs_metrics_comparison(
#       mcnemar_path=args.mcnemar,
#       perf_path=args.perf,
#       output_path=args.output,
#       dpi=args.dpi,
#       verbose=True
#   )
  
    plot_performance_complementarity(
        mcnemar_path=args.mcnemar,
        perf_path=args.perf,
        output_path='/Users/kerimulterer/prsInteractive/results/',
        multi_phenotype= [{'label':'type2Diabetes',
          'mcnemar_path':'/Users/kerimulterer/prsInteractive/results/type2Diabetes/combinedAnalysis/scores/McNemarStatsTestsAcrossPRSCalculations_refactored.csv',
          'perf_path':'/Users/kerimulterer/prsInteractive/results/type2Diabetes/combinedAnalysis/scores/model_recall_precision_improvement.csv',
          'models_to_keep':['main','epi_product','epi_summed','cardio_product','cardio_summed']},
          {'label':'celiacDisease',
                    'mcnemar_path':'/Users/kerimulterer/prsInteractive/results/celiacDisease/combinedAnalysis/scores/McNemarStatsTestsAcrossPRSCalculations_refactored.csv',
                    'perf_path':'/Users/kerimulterer/prsInteractive/results/celiacDisease/combinedAnalysis/scores/model_recall_precision_improvement.csv',
                    'models_to_keep':['main','epi_product','cardio_product']}]
            )
    
    plt.show()
    
#   Optional list of dicts, one per phenotype, each with keys:
#                               label          : str  — phenotype display name
#                               mcnemar_path   : str  — path to McNemar CSV
#                               perf_path      : str  — path to perf CSV
#                               models_to_keep : list — models to include (or None)
      