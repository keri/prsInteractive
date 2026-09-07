#!/usr/bin/env python3
# Keri Multerer 2026
# Stage 1 of the two-stage feature reduction: cheap univariate screen that
# prunes the candidate feature space before any model is fit.
#
# The screen reads pre-computed additive summary statistics (SAIGE, GRG
# `grapp assoc`, or any fastGWA-format file) rather than touching the .raw
# genotype matrices, so it costs seconds instead of the hours the per-section
# model fits take. Only variants that survive the screen are handed to
# stage 2 (sklearnSectionModelsScoreTrainV2.py).
#
# Epistasis note: a pair is retained if EITHER member survives, so that
# interaction effects whose individual marginal signal is weak are not
# discarded by an additive filter.

import argparse
import os

import pandas as pd


def load_sumstats(path):
    """fastGWA-style summary statistics; needs at least SNP and P."""
    df = pd.read_csv(path, sep=None, engine='python')
    rename = {'ID': 'SNP', 'MarkerID': 'SNP', 'p.value': 'P', 'P_BOLT_LMM': 'P'}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    missing = {'SNP', 'P'} - set(df.columns)
    if missing:
        raise ValueError(f'summary statistics missing required columns: {missing}')
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    return df.dropna(subset=['P'])


def screen_snps(sumstats, p_threshold):
    return set(sumstats.loc[sumstats['P'] < p_threshold, 'SNP'].astype(str))


def screen_snplist(input_file, kept, out_file):
    """Filter a snplist (one SNP per line) or a PLINK .bim (SNP in column 2)."""
    snps = []
    with open(input_file) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            # .bim/.map: CHR SNP CM POS A1 A2 -> the id is the second field
            snps.append(parts[1] if len(parts) >= 4 else parts[0])
    survivors = [s for s in snps if s in kept]
    with open(out_file, 'w') as f:
        f.write('\n'.join(survivors) + '\n')
    return len(snps), len(survivors)


def screen_epi_file(input_file, kept, out_file):
    """Filter a PLINK epistasis summary file, keeping pairs where either
    member survives the univariate screen."""
    df = pd.read_csv(input_file, sep=r'\s+')
    if not {'SNP', 'BEST_SNP'}.issubset(df.columns):
        raise ValueError(f'{input_file} does not look like a plink .epi summary file')
    mask = df['SNP'].astype(str).isin(kept) | df['BEST_SNP'].astype(str).isin(kept)
    out = df[mask]
    out.to_csv(out_file, sep=' ', index=False)
    return len(df), len(out)


def main():
    parser = argparse.ArgumentParser(
        description='Univariate pre-screen for the two-stage feature reduction')
    parser.add_argument('--sumstats', required=True,
                        help='fastGWA/SAIGE/GRG summary statistics file')
    parser.add_argument('--input_file', required=True,
                        help='snplist (main) or plink epi summary file (epi)')
    parser.add_argument('--data_type', required=True, choices=['main', 'epi'])
    parser.add_argument('--out_file', required=True)
    parser.add_argument('--p_threshold', type=float, default=1e-3)
    parser.add_argument('--max_features', type=int, default=None,
                        help='cap on retained features; keeps the most significant')
    args = parser.parse_args()

    sumstats = load_sumstats(args.sumstats)
    print(f'loaded {len(sumstats)} variants from {args.sumstats}')

    if args.max_features:
        ranked = sumstats.nsmallest(args.max_features, 'P')
        kept = set(ranked['SNP'].astype(str))
        effective_p = ranked['P'].max()
        print(f'keeping the top {args.max_features} variants '
              f'(effective P threshold {effective_p:.3g})')
    else:
        kept = screen_snps(sumstats, args.p_threshold)
        print(f'{len(kept)} variants pass P < {args.p_threshold:g}')

    os.makedirs(os.path.dirname(args.out_file) or '.', exist_ok=True)

    if args.data_type == 'main':
        n_in, n_out = screen_snplist(args.input_file, kept, args.out_file)
        unit = 'SNPs'
    else:
        n_in, n_out = screen_epi_file(args.input_file, kept, args.out_file)
        unit = 'pairs'

    reduction = (1 - n_out / n_in) * 100 if n_in else 0
    print(f'{unit}: {n_in} -> {n_out} ({reduction:.1f}% reduction) -> {args.out_file}')
    if n_out == 0:
        print('WARNING: nothing survived the screen; loosen --p_threshold')


if __name__ == '__main__':
    main()
