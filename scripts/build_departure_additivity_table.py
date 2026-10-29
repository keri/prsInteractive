#!/usr/bin/env python3
# Keri Multerer 2026
# Supplementary table for the departure-from-additivity analysis, combined
# across phenotypes and written as an Excel workbook.
#
# SCOPE
# Only features carrying a non-zero coefficient in the final ePRS calculations
# are included, matching what the manuscript reports. The joint-model runs test
# every decoupled candidate, but most are discarded by L1 penalisation and never
# enter a score, so reporting them would describe features the ePRS does not
# use.
#
# For the same reason CD contributes only its G+GxE environmental features:
# GxGxE was excluded as redundant in that disease (Table 2), so its terms are
# not part of any CD score. T2D retains both environmental encodings and
# contributes both.
#
# MULTIPLE TESTING
# Benjamini-Hochberg is recomputed within each reported set rather than
# inherited from the full candidate run. Restricting the features changes the
# number of tests, so a q-value carried over from the larger run would not
# correspond to the set being reported. This is why the CD environmental count
# is 10 rather than the 11 obtained when the two encodings are pooled.
#
# OUTPUT
# Extended Table S16: one workbook with a notes sheet and one sheet per
# analysis, each carrying a phenotype column.

import argparse
import os

import numpy as np
import pandas as pd

PRS = '/Users/kerimulterer/prsInteractive'
SNP_MODELS = [('epi_product', 'PRS_GxG'), ('epi_summed', 'PRS_G+G')]
# (phenotype, encoding, model column, label) -- CD omits GxGxE as redundant
GXE_SETS = [
    ('type2Diabetes', 'summed',  'cardio_summed',  'PRS_G+GxE'),
    ('type2Diabetes', 'product', 'cardio_product', 'PRS_GxGxE'),
    ('celiacDisease', 'summed',  'cardio_summed',  'PRS_G+GxE'),
]

NOTES = [
    ('Scope', 'Only features with a non-zero coefficient in the final ePRS '
              'calculations are included. The joint-model runs test every decoupled '
              'candidate, but most are dropped by L1 penalisation and never enter a '
              'score.'),
    ('CD environmental features', 'CD contributes only G+GxE. GxGxE was excluded as '
              'redundant in coeliac disease (Table 2), so its terms are not used in any '
              'CD score. T2D retains both encodings and contributes both.'),
    ('departure_from_additive', 'TRUE when q_value < 0.05, i.e. the pair or term shows '
              'significant departure from additivity within its reported set. Counts of TRUE '
              'reproduce the figures in the manuscript.'),
    ('q_value', 'Benjamini-Hochberg, recomputed within each reported set rather than '
              'inherited from the full candidate run. Restricting the features changes '
              'the number of tests, so a carried-over q would not correspond to the set '
              'reported here.'),
    ('J', 'Coefficient on the product term in logit(y) ~ C + b1*g1 + b2*g2 + J*(g1*g2). '
              'It is a PARTIAL effect: what the product adds once each locus has its own '
              'additive effect in the model. ePRS coefficients are different quantities '
              '-- total effects of the encoded feature, fitted without main effects.'),
    ('Interpretation', 'A significant J means the two loci interact. It does not mean '
              'the encoded ePRS feature is non-additive. For G+G the encoding is a sum '
              'and so additive by construction, whatever J shows.'),
    ('testable', 'Pairs or terms with enough double carriers to fit the joint model '
              '(>2,000 individuals). Excluded features have no J and are not counted in '
              'any denominator.'),
    ('n_double_carriers / n_carriers', 'Individuals carrying both constituent alleles, '
              'which determines whether a term is testable and how precisely J is '
              'estimated.'),
    ('HLA linkage', 'In CD, DQB1*02:01, DRB1*03:01 and DRB3*01:01 lie on the same DR3-B8 '
              'haplotype. Rows involving them are one signal observed through linked '
              'alleles, not independent findings.'),
]


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    q[order] = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def retained(pheno, model):
    fs = pd.read_csv(os.path.join(
        PRS, 'results', pheno, 'combinedAnalysis', 'scores',
        'featureScoresReducedFinalModel.combined.filtered.csv'))
    x = fs[(fs.model == model) & (fs.coefs != 0) & (fs.feature != '(Intercept)')]
    return dict(zip(x.feature, x.coefs))


def snp_pairs(results_root):
    rows = []
    for pheno, label in [('type2Diabetes', 'T2D'), ('celiacDisease', 'CD')]:
        ar = pd.read_csv(os.path.join(results_root, pheno, 'gxg_joint_dc2000',
                                      'gxg_joint_model_results.csv'))
        ar['key'] = [frozenset((a, b)) for a, b in zip(ar.feature1, ar.feature2)]
        for model, name in SNP_MODELS:
            ret = retained(pheno, model)
            keys = {frozenset(t.strip() for t in f.split(',')): f
                    for f in ret if f.count(',') == 1}
            sub = ar[ar.key.isin(keys) & ar.p_J.notna()].copy()
            if not len(sub):
                continue
            sub['q_value'] = bh(sub.p_J.values)
            sub['phenotype'] = label
            sub['eprs_model'] = name
            sub['feature'] = [keys[k] for k in sub.key]
            sub['eprs_coefficient'] = [ret[f] for f in sub.feature]
            rows.append(sub.drop(columns=['key']))
            print(f'  {label:<4}{name:<10} testable {len(sub):>5}  '
                  f'q<0.05 {int((sub.q_value < 0.05).sum()):>5}')
    return pd.concat(rows, ignore_index=True)


def gxe_terms(results_root):
    env = set(pd.read_csv(os.path.join(PRS, 'results', 'participant_environment.csv'),
                          nrows=0).columns) - {'Participant ID'}
    frames = {}
    for pheno, enc, model, name in GXE_SETS:
        label = 'T2D' if pheno == 'type2Diabetes' else 'CD'
        ar = pd.read_csv(os.path.join(results_root, pheno, f'gxe_joint_{enc}',
                                      'gxe_joint_model_results.csv'))
        ar['key'] = [frozenset([r['env']] + [v for v in (r.get('feature1'), r.get('feature2'))
                                             if isinstance(v, str)])
                     for _, r in ar.iterrows()]
        ret = retained(pheno, model)
        keys = {frozenset(t.strip() for t in f.split(',')): f for f in ret
                if any(t.strip() in env for t in f.split(','))}
        sub = ar[ar.key.isin(keys) & ar.p_J.notna()].copy()
        sub['phenotype'] = label
        sub['eprs_model'] = name
        sub['feature'] = [keys[k] for k in sub.key]
        sub['eprs_coefficient'] = [ret[f] for f in sub.feature]
        frames.setdefault(label, []).append(sub.drop(columns=['key']))

    out = []
    # BH within each phenotype's reported set: T2D pools both encodings (41
    # features), CD has only G+GxE (102).
    for label, fr in frames.items():
        d = pd.concat(fr, ignore_index=True)
        d['q_value'] = bh(d.p_J.values)
        print(f'  {label:<4}{"GxE":<10} testable {len(d):>5}  '
              f'q<0.05 {int((d.q_value < 0.05).sum()):>5}')
        out.append(d)
    return pd.concat(out, ignore_index=True)


# q_J and non_additive come from the original run, where BH spanned every
# decoupled candidate. Restricting to non-zero features changes the number of
# tests, so those values do not apply to the set reported here and are dropped
# rather than left beside the recomputed q_value for a reader to confuse.
STALE = ['q_J', 'non_additive']


Q_THRESHOLD = 0.05


def tidy(d, first):
    d = d.drop(columns=[c for c in STALE if c in d.columns])
    # Derived from the recomputed q_value, so it matches the counts reported in
    # the manuscript rather than the dropped flag from the full candidate run.
    d['departure_from_additive'] = d['q_value'] < Q_THRESHOLD
    cols = [c for c in first if c in d.columns] + \
           [c for c in d.columns if c not in first]
    return d[cols].sort_values(['phenotype', 'eprs_model', 'q_value'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results_root', default='results')
    ap.add_argument('--out', required=True, help='output .xlsx path')
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    print('SNP-pair features (G+G and GxG)')
    pairs = tidy(snp_pairs(args.results_root),
                 ['phenotype', 'eprs_model', 'feature', 'feature1', 'feature2',
                  'eprs_coefficient', 'n', 'n_double_carriers', 'J', 'se_J',
                  'p_J', 'q_value', 'departure_from_additive'])
    print('\nEnvironmental interaction features (GxE)')
    gxe = tidy(gxe_terms(args.results_root),
               ['phenotype', 'eprs_model', 'feature', 'env', 'feature1', 'feature2',
                'order', 'eprs_coefficient', 'n', 'n_carriers', 'J', 'se_J',
                'p_J', 'q_value', 'departure_from_additive'])

    with pd.ExcelWriter(args.out, engine='openpyxl') as w:
        pd.DataFrame(NOTES, columns=['Item', 'Definition']).to_excel(
            w, sheet_name='Notes', index=False)
        pairs.to_excel(w, sheet_name='SNP_pair_features', index=False)
        gxe.to_excel(w, sheet_name='GxE_features', index=False)
        for sheet, width in (('Notes', (28, 110)), ('SNP_pair_features', (34, 16)),
                             ('GxE_features', (40, 16))):
            ws = w.sheets[sheet]
            ws.column_dimensions['A'].width = width[0]
            ws.column_dimensions['B'].width = width[1]
            ws.freeze_panes = 'A2'
    print(f'\n-> {args.out}')
    print(f'   SNP_pair_features {len(pairs)} rows, GxE_features {len(gxe)} rows')


if __name__ == '__main__':
    main()
