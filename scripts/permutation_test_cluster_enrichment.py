#!/usr/bin/env python3
# Keri Multerer 2026
# Permutation test for the differential correspondence between ePRS feature
# sets and T2D functional clusters / disease subtypes.
#
# Figure 4 shows concordance-weighted feature counts per model, per functional
# cluster and per subtype. Those counts are descriptive: a model contributing
# many features to one cluster may simply have more features overall, or the
# cluster may be large. This tests whether the concentration is greater than
# expected if features were exchangeable between models.
#
# NULL MODEL
# Each feature carries its own cluster/subtype annotations and concordance
# values, which are properties of the feature rather than of the model that
# selected it (verified: concordance never varies by cohort for a given
# feature/cluster/subtype). The null therefore permutes the model label across
# features, holding each feature's annotation profile intact and preserving the
# number of features per model. This breaks only the model-feature association,
# which is the thing being tested.
#
# STATISTIC
# Share = the number of a model's distinct features assigned to that cluster
# (or subtype), divided by the model's total assigned features. Shares rather
# than raw counts, so models with more features are comparable. Counts rather
# than summed concordance, because concordance is signed and one model's total
# is negative, which makes concordance shares uninterpretable.
#
# Two-sided empirical p-values with a +1 correction, then Benjamini-Hochberg
# across all tested cells.
#
# NOTE ON POWER: models carry 11-27 features each, so cells are small and the
# minimum attainable p-value is bounded by --n_perm. Absence of significance
# should be read as insufficient evidence, not as evidence of absence.

import argparse
import os

import numpy as np
import pandas as pd


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.clip(adj, 0, 1)
    return q


def weighted_shares(assign, group_col):
    """Share of each model's features that fall in each group.

    The statistic counts distinct features rather than summing concordance.
    Concordance is signed (concordant 1, discordant -1, neutral 0), so a
    model whose discordant features outweigh its concordant ones has a
    negative total: cardio_product sums to -8. Dividing a signed subtotal by a
    signed total produces uninterpretable shares, including positive shares
    from negative counts. Feature counts are strictly non-negative and answer
    the question directly -- does this model contribute a disproportionate
    share of its features to this cluster. Net concordance is carried through
    as a descriptive column.
    """
    g = (assign.groupby(['cohort', group_col])
         .agg(n_features=('feature', 'nunique'),
              net_concordance=('concordance', 'sum'))
         .reset_index())
    tot = (assign.groupby('cohort').feature.nunique()
           .rename('cohort_total').reset_index())
    g = g.merge(tot, on='cohort')
    g['share'] = g.n_features / g.cohort_total
    return g


def run_permutation(assign, group_col, n_perm, seed):
    obs = weighted_shares(assign, group_col)
    key = ['cohort', group_col]

    # one row per (cohort, feature); features keep their annotation rows
    pairs = assign[['cohort', 'feature']].drop_duplicates().reset_index(drop=True)
    feature_rows = assign.drop(columns=['cohort'])

    rng = np.random.default_rng(seed)
    labels = pairs.cohort.to_numpy()
    ge = np.zeros(len(obs)); le = np.zeros(len(obs))
    obs_idx = {tuple(r): i for i, r in enumerate(obs[key].to_numpy())}

    for _ in range(n_perm):
        permuted = pairs.copy()
        permuted['cohort'] = rng.permutation(labels)
        sim = permuted.merge(feature_rows, on='feature', how='left')
        s = weighted_shares(sim, group_col)
        vals = np.zeros(len(obs))
        for r, v in zip(s[key].to_numpy(), s.share.to_numpy()):
            i = obs_idx.get(tuple(r))
            if i is not None:
                vals[i] = v
        ge += vals >= obs.share.to_numpy()
        le += vals <= obs.share.to_numpy()

    obs['p_greater'] = (ge + 1) / (n_perm + 1)
    obs['p_less'] = (le + 1) / (n_perm + 1)
    obs['p_two_sided'] = np.minimum(1.0, 2 * np.minimum(obs.p_greater, obs.p_less))
    obs['q_value'] = benjamini_hochberg(obs.p_two_sided.to_numpy())
    obs['direction'] = np.where(obs.p_greater < obs.p_less, 'over', 'under')
    obs['n_perm'] = n_perm
    return obs.sort_values('p_two_sided')


def main():
    ap = argparse.ArgumentParser(
        description='Permutation test for ePRS feature set correspondence with '
                    'functional clusters and T2D subtypes')
    ap.add_argument('--assignments_file', required=True,
                    help='functionalEnrichment/subtype_feature_assignments.csv')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--n_perm', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--q_threshold', type=float, default=0.05)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    assign = pd.read_csv(args.assignments_file)
    required = {'feature', 'cohort', 'cluster_key', 'subtype', 'concordance'}
    missing = required - set(assign.columns)
    if missing:
        raise SystemExit(f'assignments file missing columns: {missing}')
    assign = assign.dropna(subset=['cohort', 'feature', 'concordance'])

    n_feat = assign.groupby('cohort').feature.nunique()
    print(f'{len(assign)} annotation rows, {assign.feature.nunique()} unique features, '
          f'{assign.cohort.nunique()} models')
    print('features per model:')
    print(n_feat.to_string())
    print(f'\nrunning {args.n_perm} permutations (seed {args.seed}) ...\n')

    for group_col, label in [('cluster_key', 'functional cluster'), ('subtype', 'T2D subtype')]:
        res = run_permutation(assign, group_col, args.n_perm, args.seed)
        out = os.path.join(args.out_dir, f'permutation_enrichment_{group_col}.csv')
        res.to_csv(out, index=False)
        sig = res[res.q_value < args.q_threshold]
        print(f'=== {label}: {len(res)} model-by-{group_col} cells tested, '
              f'{len(sig)} with q < {args.q_threshold} ===')
        show = sig if len(sig) else res.head(5)
        cols = ['cohort', group_col, 'n_features', 'cohort_total', 'share',
                'net_concordance', 'p_two_sided', 'q_value', 'direction']
        print(show[cols].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
        if not len(sig):
            print(f'(none significant; showing the five smallest p-values. '
                  f'Minimum attainable p is {1/(args.n_perm+1):.2g})')
        print(f'-> {out}\n')


if __name__ == '__main__':
    main()
