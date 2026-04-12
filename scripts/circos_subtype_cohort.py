#!/usr/bin/env python3
"""
circos_subtype_cohort.py
========================
One circos-style bar plot per cohort from subtype_feature_assignments.csv.

Layout:
  • Genomic (loci): only OR > 1 shown; bars extend OUTWARD; labelled with gene name
  • Clinical:       effect_size_r > 0 → bar OUTWARD; < 0 → bar INWARD (below midline)
  • Bar colour:     red = OR>1 or r>0 (risk ↑),  blue = r<0 (protective)
  • Baseline ring:  segmented green arcs, one per cluster with gap between clusters
  • Cluster label:  inside the ring at cluster midpoint, black, tangential
  • Feature label:  at bar tip, 90° to bar (tangential), black

Usage
-----
python circos_subtype_cohort.py \\
    --assignments results/.../subtype_feature_assignments.csv \\
    --out_dir     results/.../enrichment/plots
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── colours ───────────────────────────────────────────────────────────────────
BAR_POS      = '#C0392B'   # red  – genomic OR > 1
BAR_NEG      = '#2471A3'   # blue – genomic (unused after filtering, kept for fallback)
CLINICAL_CLR = '#C0B52B'   # yellow-green – clinical bars (solid r>0, hatched r<0)

SUBTYPE_COLORS = {
    'SAID': '#4472C4',
    'SIDD': '#ED7D31',
    'SIRD': '#9B59B6',
    'MOD':  '#2ECC71',
    'MARD': '#F1C40F',
}

CLUSTER_ORDER = [
    'beta_cell_1', 'beta_cell_2', 'proinsulin', 'hyper_insulin',
    'res_glycaemic', 'body_fat', 'obesity', 'metabolic_syndrome',
    'lipodystrophy', 'liver_lipid',
]

CLUSTER_LABELS = {
    'beta_cell_1':        'Beta Cell 1',
    'beta_cell_2':        'Beta Cell 2',
    'proinsulin':         'Proinsulin',
    'hyper_insulin':      'Hyper Insulin',
    'res_glycaemic':      'Res. Glycaemic',
    'body_fat':           'Body Fat',
    'obesity':            'Obesity',
    'metabolic_syndrome': 'Metab. Syndrome',
    'lipodystrophy':      'Lipodystrophy',
    'liver_lipid':        'Liver / Lipid',
}

# 3-letter abbreviations shown on the arc; full names go in the legend
CLUSTER_ABBREV = {
    'beta_cell_1':        'BC1',
    'beta_cell_2':        'BC2',
    'proinsulin':         'PRO',
    'hyper_insulin':      'HIN',
    'res_glycaemic':      'RGY',
    'body_fat':           'BFT',
    'obesity':            'OBS',
    'metabolic_syndrome': 'MET',
    'lipodystrophy':      'LIP',
    'liver_lipid':        'LVR',
}

CLUSTER_COLORS = {
    'beta_cell_1':        '#27AE60',
    'beta_cell_2':        '#1ABC9C',
    'proinsulin':         '#E67E22',
    'hyper_insulin':      '#D35400',
    'res_glycaemic':      '#C0392B',
    'body_fat':           '#3498DB',
    'obesity':            '#2980B9',
    'metabolic_syndrome': '#16A085',
    'lipodystrophy':      '#8E44AD',
    'liver_lipid':        '#5B2C6F',
}

# ── helpers ───────────────────────────────────────────────────────────────────

def load_and_prepare(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[~df['cluster_key'].fillna('').str.startswith('direct__')].copy()


def _is_real_gene(g: str) -> bool:
    """Return True if the symbol looks like a real named gene (not LOC/LINC/pseudogene)."""
    if not g or g.lower() == 'nan':
        return False
    u = g.upper()
    if re.match(r'^LOC\d',  u): return False   # e.g. LOC105378657
    if re.match(r'^LINC\d', u): return False   # e.g. LINC01234
    if re.match(r'^MIR\d',  u): return False   # microRNA
    if re.match(r'^ENSG\d', u): return False   # raw Ensembl ID
    if u.endswith('-AS1') or u.endswith('-AS2'): return False  # antisense
    return True


def _parse_genes(val) -> list:
    """Split a semicolon/comma-separated gene string into individual symbols."""
    if not pd.notna(val) or str(val).strip().lower() in ('', 'nan'):
        return []
    return [g.strip() for g in re.split(r'[;,\s]+', str(val).strip())
            if g.strip() and g.strip().lower() != 'nan']


def _best_gene(snp1, snp2) -> str:
    """
    Pick the most biologically meaningful gene symbol.

    Priority:
      1. Real genes from snp2 (the nominated / causal gene column)
      2. Real genes from snp1 (the flanking gene column)
      3. Any gene from snp2 (even LOC)
      4. Any gene from snp1
    """
    snp2_real = [g for g in _parse_genes(snp2) if _is_real_gene(g)]
    snp1_real = [g for g in _parse_genes(snp1) if _is_real_gene(g)]

    # Prefer snp2 single real gene — most likely the causal/nominated locus gene
    if len(snp2_real) == 1:
        return snp2_real[0]
    # snp2 has multiple real genes — take first
    if snp2_real:
        return snp2_real[0]
    # Fall back to snp1 real genes
    if snp1_real:
        return snp1_real[0]
    # Last resort — any token from snp2 then snp1
    for val in (snp2, snp1):
        g = _parse_genes(val)
        if g:
            return g[0]
    return ''


def _is_rsid(feat: str) -> bool:
    """Return True if the feature looks like an rs-ID (genetic locus)."""
    return bool(re.match(r'^rs\d+$', str(feat).strip(), re.IGNORECASE))


def _gene_label(row) -> str:
    """Return best gene symbol for genomic/rsID features; cleaned name for clinical traits."""
    feat = str(row.get('feature', ''))
    # Always use gene lookup for rs-IDs regardless of feature_source label
    if _is_rsid(feat) or str(row.get('feature_source', '')).strip() != 'clinical':
        best = _best_gene(row.get('snp1_gene_names'), row.get('snp2_gene_names'))
        return best if best else feat
    # True clinical trait: clean up the name
    feat = re.sub(r'[_\s]+', ' ', feat).strip()
    feat = re.sub(r'\s*\(.*\)', '', feat)
    return feat.title()


def build_bar_frame(df: pd.DataFrame, cohort: str,
                    include_protective: bool = True) -> pd.DataFrame:
    sub = df[df['cohort'] == cohort].copy()
    if sub.empty:
        return pd.DataFrame()

    gcols = ['feature', 'cluster_key', 'display_name', 'feature_source']
    for c in ('snp1_gene_names', 'snp2_gene_names', 'training_data_used'):
        if c in sub.columns:
            gcols.append(c)

    agg = (sub.groupby(gcols, sort=False, dropna=False)
              .agg(odds_ratio=('odds_ratio', 'mean'),
                   effect_size_r=('effect_size_r', 'mean'))
              .reset_index())

    def _signed(row):
        if pd.notna(row.get('effect_size_r')) and row.get('effect_size_r') != 0:
            return float(row['effect_size_r'])
        if pd.notna(row.get('odds_ratio')):
            return float(np.log2(row['odds_ratio']))
        return 0.0

    agg['bar_value'] = agg.apply(_signed, axis=1)

    # rs-IDs are genomic loci regardless of feature_source label (some are tagged
    # 'clinical' via the Mahajan bridge but are still genetic variants).
    is_rsid = agg['feature'].str.match(r'^rs\d+$', case=False, na=False)
    is_genomic = is_rsid | (agg['feature_source'] == 'genomic_shap')

    # LowControls features represent the protective cohort signature:
    # OR < 1 → bar_value < 0 → shown as inward (blue) bars.
    # HighCases / ExclusiveHighCases: only OR > 1 (risk direction).
    is_lowcontrols = agg.get('training_data_used', pd.Series('', index=agg.index)) == 'LowControls'

    keep_genomic_risk = is_genomic & ~is_lowcontrols & (agg['bar_value'] > 0)

    if include_protective:
        keep_genomic_protective = is_genomic & is_lowcontrols & (agg['bar_value'] < 0)
        # Clinical traits: both directions
        keep_clinical = (~is_rsid) & (agg['feature_source'] == 'clinical')
    else:
        # Risk-only plot: drop all LowControls and negative-r clinical features
        keep_genomic_protective = pd.Series(False, index=agg.index)
        keep_clinical = (~is_rsid) & (agg['feature_source'] == 'clinical') & (agg['bar_value'] > 0)

    agg = pd.concat([
        agg[keep_genomic_risk],
        agg[keep_genomic_protective],
        agg[keep_clinical],
    ]).reset_index(drop=True)

    agg['bar_height'] = agg['bar_value'].abs()
    agg['is_inward']  = agg['bar_value'] < 0

    # Colour / hatch scheme:
    #   Genomic risk     (OR > 1, outward) : solid red      BAR_POS
    #   Genomic protect  (OR < 1, inward)  : hashed red     BAR_POS + hatch '//'
    #   Clinical risk    (r > 0, outward)  : solid clinical  CLINICAL_CLR
    #   Clinical protect (r < 0, inward)   : hashed clinical CLINICAL_CLR + hatch '//'
    is_rsid_final     = agg['feature'].str.match(r'^rs\d+$', case=False, na=False)
    is_hla_final      = agg['feature'].str.match(r'^HLA[-_*]', case=False, na=False)
    is_genomic_final  = is_rsid_final | is_hla_final | (agg['feature_source'] == 'genomic_shap')
    is_clinical_final = ~is_genomic_final

    agg['bar_color'] = np.where(is_clinical_final, CLINICAL_CLR, BAR_POS)
    agg['bar_hatch'] = np.where(agg['bar_value'] < 0, '//', '')
    agg['label']     = agg.apply(_gene_label, axis=1)
    return agg


def build_subtype_map(df: pd.DataFrame, cohort: str) -> dict:
    """
    Returns dict keyed by (feature, cluster_key) → {'subtypes': [...], 'concordance': float}

    concordance is the mean concordance value for this feature×cluster pair
    (+1 concordant, -1 discordant, 0 neutral).  Falls back gracefully when the
    column is absent (old CSV files without concordance).
    """
    sub = df[df['cohort'] == cohort]
    has_conc = 'concordance' in sub.columns
    out = {}
    for (feat, ck), grp in sub.groupby(['feature', 'cluster_key']):
        subtypes = sorted(grp['subtype'].dropna().unique().tolist())
        conc = float(grp['concordance'].mean()) if has_conc else 0.0
        out[(feat, ck)] = {'subtypes': subtypes, 'concordance': conc}
    return out


def ordered_bars(agg: pd.DataFrame) -> list:
    present = [c for c in CLUSTER_ORDER if c in agg['cluster_key'].values]
    extra   = [c for c in agg['cluster_key'].unique() if c not in present]
    bars = []
    for ck in present + sorted(extra):
        chunk = agg[agg['cluster_key'] == ck].sort_values('bar_height', ascending=False)
        # Deduplicate by gene label within each cluster — keep highest bar per label
        seen_labels: set = set()
        for rec in chunk.to_dict('records'):
            lbl = rec.get('label', rec.get('feature', ''))
            if lbl not in seen_labels:
                seen_labels.add(lbl)
                bars.append(rec)
    return bars


def assign_angles(bars: list, gap_deg: float) -> tuple:
    """CW-from-top degree angle per bar centre, with gaps between clusters."""
    keys = list(dict.fromkeys(b['cluster_key'] for b in bars))
    dpb  = (360.0 - gap_deg * len(keys)) / max(len(bars), 1)
    angles, cursor, prev = [], gap_deg / 2, None
    for b in bars:
        ck = b['cluster_key']
        if prev is not None and ck != prev:
            cursor += gap_deg
        angles.append(cursor + dpb / 2)
        cursor += dpb
        prev = ck
    return np.array(angles), dpb


def _radial_rot(a_deg: float) -> float:
    """
    Rotation (°) so text runs radially — parallel to the bar, reading outward.
    Flip only when the raw rotation would produce upside-down characters,
    i.e. when a_deg > 180 (left half of circle).
    """
    rot = 90.0 - a_deg
    if (a_deg % 360) > 180:
        rot += 180.0
    return rot


# ── main draw ─────────────────────────────────────────────────────────────────

def draw_circos(bars: list, subtype_map: dict, cohort: str, out_path: str,
                show_labels: bool = True, include_protective: bool = True):
    if not bars:
        print(f'  [{cohort}] no bars – skipping')
        return

    n       = len(bars)
    # Spread sparse cohorts around full circle; dense cohorts use smaller gaps
    gap_deg = max(5.0, min(22.0, 500.0 / max(n, 1)))
    figsize = 22 if n > 80 else (17 if n > 30 else 14)

    # ── radii ──
    R_BASE   = 1.0
    MAX_OUT  = 0.44    # max outward bar height
    MAX_IN   = 0.22    # max inward (clinical r<0) height
    R_DOT_0  = R_BASE + MAX_OUT + 0.07
    DOT_STP  = 0.068
    N_DOT    = 5
    R_CLBL   = R_BASE - 0.16    # cluster labels inside ring

    angles_deg, dpb = assign_angles(bars, gap_deg)
    angles_rad      = np.radians(90.0 - angles_deg)

    # Normalise heights (outward and inward separately)
    h_raw    = np.array([b['bar_height'] for b in bars])
    out_mask = np.array([not b.get('is_inward', False) for b in bars])
    in_mask  = ~out_mask

    vmax_out = h_raw[out_mask].max() if out_mask.any() else 1.0
    vmax_in  = h_raw[in_mask].max()  if in_mask.any()  else 1.0
    if vmax_out == 0: vmax_out = 1.0
    if vmax_in  == 0: vmax_in  = 1.0

    heights  = np.where(out_mask,
                        h_raw / vmax_out * MAX_OUT,
                        h_raw / vmax_in  * MAX_IN)

    # ── figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(figsize, figsize), facecolor='white')
    ax  = fig.add_subplot(111, projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_axis_off()
    ax.set_rlim(0, R_DOT_0 + N_DOT * DOT_STP + 0.65)

    # ── cluster spans ─────────────────────────────────────────────────────────
    spans: dict = {}
    for b, a in zip(bars, angles_deg):
        ck = b['cluster_key']
        if ck not in spans:
            spans[ck] = {'s': a, 'e': a}
        else:
            spans[ck]['e'] = a

    # ── segmented black baseline arcs (just short of edge bar centres) ───────
    arc_pad = dpb * 0.08   # tiny gap so arc ends just before edge node centres
    for ck, sp in spans.items():
        half = dpb / 2
        a0   = np.radians(90 - (sp['s'] - half + arc_pad))
        a1   = np.radians(90 - (sp['e'] + half - arc_pad))
        npt  = max(int(abs(a0 - a1) * 200) + 2, 8)
        t_arc = np.linspace(a0, a1, npt)
        ax.plot(t_arc, np.full(npt, R_BASE),
                color='black', lw=3.0, zorder=3, solid_capstyle='butt')

    # ── bars ──────────────────────────────────────────────────────────────────
    bw = np.radians(dpb * 0.68)
    for b, th, h in zip(bars, angles_rad, heights):
        hatch = b.get('bar_hatch', '')
        if b.get('is_inward', False):
            ax.bar(th, h, width=bw, bottom=R_BASE - h,
                   color=b['bar_color'], hatch=hatch, alpha=0.85, zorder=4, linewidth=0)
        else:
            ax.bar(th, h, width=bw, bottom=R_BASE,
                   color=b['bar_color'], hatch=hatch, alpha=0.85, zorder=4, linewidth=0)

    # ── cluster 3-letter abbreviations: grey, just inside the ring, no overlap ─
    # Placed at R_BASE - small offset so they sit inside the green arc away from bars
    cl_fs = max(6.5, min(9.5, 13.0 - n // 25))
    for ck, sp in spans.items():
        abbrev  = CLUSTER_ABBREV.get(ck, ck[:3].upper())
        mid_cw  = (sp['s'] + sp['e']) / 2
        mid_pol = np.radians(90 - mid_cw)
        rot     = _radial_rot(mid_cw)
        r_clbl  = R_BASE - 0.13   # inside the ring, well below any outward bar
        ax.text(mid_pol, r_clbl, abbrev,
                ha='center', va='center',
                fontsize=cl_fs, fontweight='normal',
                color='#AAAAAA', rotation=rot, rotation_mode='anchor',
                zorder=6)

    # ── subtype dots: concordance-weighted, one set per cluster at midpoint ──
    dot_sz = max(14, min(50, 1300 / max(len(spans), 1)))
    for ck, sp in spans.items():
        mid_cw  = (sp['s'] + sp['e']) / 2
        mid_pol = np.radians(90 - mid_cw)
        cluster_bars = [b for b in bars if b['cluster_key'] == ck]

        # Accumulate concordance-weighted score per subtype.
        # concordance ∈ {+1, 0, -1}; neutral (0) features contribute a small
        # positive weight (0.1) so they are not silenced but can be overridden
        # by explicitly discordant evidence.
        subtype_scores: dict = {}
        for b in cluster_bars:
            info = subtype_map.get((b['feature'], ck), {})
            if isinstance(info, dict):
                subtypes = info.get('subtypes', [])
                conc     = info.get('concordance', 0.0)
            else:
                subtypes = info          # legacy list fallback
                conc     = 0.0
            weight = conc if conc != 0.0 else 0.1
            for st in subtypes:
                subtype_scores[st] = subtype_scores.get(st, 0.0) + weight

        # Only show subtypes with net positive concordance score
        positive_subtypes = {st: sc for st, sc in subtype_scores.items() if sc > 0}
        if not positive_subtypes and subtype_scores:
            # All evidence discordant — show with low alpha as a warning marker
            positive_subtypes = subtype_scores

        max_sc = max(positive_subtypes.values()) if positive_subtypes else 1.0
        sorted_subtypes = sorted(positive_subtypes, key=lambda s: positive_subtypes[s],
                                 reverse=True)

        for j, st in enumerate(sorted_subtypes[:N_DOT]):
            sc    = positive_subtypes[st]
            alpha = float(np.clip(0.35 + 0.65 * (sc / max_sc), 0.2, 1.0))
            ax.scatter(mid_pol, R_DOT_0 + j * DOT_STP,
                       s=dot_sz,
                       color=SUBTYPE_COLORS.get(st, '#AAAAAA'),
                       alpha=alpha,
                       zorder=6, linewidths=0.5, edgecolors='white')

    # ── feature labels (skipped for clean / no-label version) ────────────────
    if show_labels:
        fs      = max(5.5, min(8.5, 320 / n))
        lbl_gap = 0.03

        for b, a_deg_val, th, h in zip(bars, angles_deg, angles_rad, heights):
            lbl = b.get('label', b.get('feature', ''))
            if len(lbl) > 20:
                lbl = lbl[:18] + '…'

            rot = _radial_rot(a_deg_val)
            ha  = 'right' if (a_deg_val % 360) > 180 else 'left'

            if b.get('is_inward', False):
                r_lbl = R_BASE - h - lbl_gap
                ha    = 'left' if (a_deg_val % 360) > 180 else 'right'
            else:
                r_lbl = R_BASE + h + lbl_gap

            ax.text(th, r_lbl, lbl,
                    ha=ha, va='center',
                    fontsize=fs, color='black',
                    rotation=rot, rotation_mode='anchor',
                    zorder=6)

    # ── centre metric note ────────────────────────────────────────────────────
    has_er = any(pd.notna(b.get('effect_size_r')) and b.get('effect_size_r', 0) != 0
                 for b in bars)
    note = (f'max |r| = {vmax_in:.3f}' if has_er
            else f'max log₂OR = {vmax_out:.3f}')
    ax.text(0, 0, note, ha='center', va='center',
            fontsize=8, color='#BBBBBB', zorder=5)

    # ── legend: colour key ────────────────────────────────────────────────────
    colour_handles = [
        mpatches.Patch(color=BAR_POS,                   label='Genomic OR > 1  (outward)'),
        mpatches.Patch(color=BAR_POS,      hatch='//',  label='Genomic OR < 1  (inward, hashed)'),
        mpatches.Patch(color=CLINICAL_CLR,              label='Clinical r > 0  (outward)'),
        mpatches.Patch(color=CLINICAL_CLR, hatch='//',  label='Clinical r < 0  (inward, hashed)'),
    ] + [mpatches.Patch(color=c, label=st) for st, c in SUBTYPE_COLORS.items()]

    leg1 = ax.legend(handles=colour_handles,
                     loc='upper right', bbox_to_anchor=(1.42, 1.10),
                     fontsize=9.5, framealpha=0.96, edgecolor='#cccccc',
                     title='Colour key', title_fontsize=9.5)
    ax.add_artist(leg1)

    # ── legend: cluster abbreviations ────────────────────────────────────────
    present_clusters = list(dict.fromkeys(b['cluster_key'] for b in bars))
    abbrev_handles = [
        mpatches.Patch(color=CLUSTER_COLORS.get(ck, '#888888'),
                       label=f"{CLUSTER_ABBREV.get(ck, ck[:3].upper())}  {CLUSTER_LABELS.get(ck, ck)}")
        for ck in present_clusters
    ]
    ax.legend(handles=abbrev_handles,
              loc='upper right', bbox_to_anchor=(1.42, 0.62),
              fontsize=9.5, framealpha=0.96, edgecolor='#cccccc',
              title='Cluster key', title_fontsize=9.5)

    title_cohort = cohort.replace('_', ' ').title()
    sig_tag = '' if include_protective else '  |  risk only'
    ax.set_title(f'Cluster–Subtype Circos  |  {title_cohort}  ({n} features){sig_tag}',
                 fontsize=14, fontweight='bold', pad=35)

    fig.savefig(out_path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [{cohort}]  {n} bars  →  {out_path}')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',
                   default='results/type2Diabetes/subtypeResearch/enrichment/'
                           'subtype_feature_assignments.csv')
    p.add_argument('--out_dir',
                   default='results/type2Diabetes/subtypeResearch/enrichment/plots')
    p.add_argument('--cohorts', nargs='+')
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f'Loading {args.assignments}')
    df = load_and_prepare(args.assignments)

    cohorts = args.cohorts or sorted(df['cohort'].dropna().unique())
    print(f'Cohorts: {cohorts}\n')

    for cohort in cohorts:
        subtype_map = build_subtype_map(df, cohort)

        for incl_prot, suffix in [(True, ''), (False, '_risk')]:
            bars = ordered_bars(build_bar_frame(df, cohort,
                                                include_protective=incl_prot))
            # labelled (Inkscape reference)
            draw_circos(bars, subtype_map, cohort,
                        os.path.join(args.out_dir,
                                     f'circos_{cohort}{suffix}_labelled.png'),
                        show_labels=True,
                        include_protective=incl_prot)
            # clean (no labels)
            draw_circos(bars, subtype_map, cohort,
                        os.path.join(args.out_dir,
                                     f'circos_{cohort}{suffix}_clean.png'),
                        show_labels=False,
                        include_protective=incl_prot)

    print('\nDone.')


if __name__ == '__main__':
    main()
