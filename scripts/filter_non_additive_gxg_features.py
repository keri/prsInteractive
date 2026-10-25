#!/usr/bin/env python3
"""
filter_non_additive_gxg_features.py

Departure-from-additivity test for retained GxG pairs — the G-by-G analogue of
prsInteractive's filter_non_additive_gen_env_features.py.

WHY THIS IS NEEDED
------------------
The ePRS pipeline runs feature selection separately on the G, G+G and GxG
feature sets, so the constituent SNPs of a retained GxG pair are almost never
also retained as main effects. Checking the decoupled coefficients for
type2Diabetes: of 2,343 distinct pair members only 3 appear as main effects,
and in 0 of 1,175 pairs are both members present.

Two consequences:

  1. A filter that compares a pair's coefficient against its constituents'
     coefficients (as the GxE filter does) would compare against absent values
     for essentially every pair, and pass everything.

  2. The pair coefficient itself is a TOTAL effect, not an interaction: with no
     main effects in the model it absorbs the marginal effects of both SNPs
     along with their interaction.

This script fits the constituents explicitly, one model per pair:

    logit(y) ~ covariates + b1*g1 + b2*g2 + J*(g1 x g2)

J is then the departure from additivity by construction, with a standard error,
a Wald p-value, and a Benjamini-Hochberg q-value across the tested pairs. J is
also the quantity the OPRS formula calls J_ij.

A minimum double-carrier requirement is applied because the product encoding is
non-zero only for individuals carrying risk alleles at both loci. In the
type2Diabetes pairs the median pair has ~14,700 double carriers, but 32 pairs
(1.5%) have fewer than 100 — and those show a median |coef| of 1.21 against
0.020 for pairs with >=1000, i.e. a ~61x inflation driven purely by low counts.

OUTPUTS (all written under --out_dir; nothing in prsInteractive is modified)
    gxg_joint_model_results.csv   every pair tested, with b1/b2/J, SEs, p, q,
                                  double-carrier counts, MAFs, convergence flag
    gxg_pairs.non_additive.csv    pairs passing the filter (the "real" GxG set)
    gxg_pairs.excluded.csv        pairs failing, with a reason column
    run_summary.txt               parameters and counts for the methods section

Example
-------
    python filter_non_additive_gxg_features.py \
        --pairs_file  results/type2Diabetes/decoupled/gxg_pairs.csv \
        --genotype_file /Users/kerimulterer/prsInteractive/results/type2Diabetes/trainingCombined.raw \
        --covar_file  /Users/kerimulterer/prsInteractive/results/covar.csv \
        --withdrawal_path /Users/kerimulterer/prsInteractive/results/type2Diabetes/withdrawals.csv \
        --out_dir     results/type2Diabetes/gxg_joint
"""

import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter(action="ignore")

COVARIATES = ["age", "SEX"] + [f"PC{i}" for i in range(1, 11)]


# ---------------------------------------------------------------- data loading

def load_pairs(path):
    """Accepts gxg_pairs.csv (feature1,feature2,...) or a per-model
    GxG_features.csv (feature1,feature2,coef)."""
    df = pd.read_csv(path)
    if not {"feature1", "feature2"}.issubset(df.columns):
        raise SystemExit(f"{path} needs feature1 and feature2 columns")
    df = df[df.feature1.notna() & df.feature2.notna()]
    df = df[(df.feature1.astype(str) != "") & (df.feature2.astype(str) != "")]
    df["feature1"] = df.feature1.astype(str)
    df["feature2"] = df.feature2.astype(str)
    # canonical ordering so a pair is not tested twice under two orientations
    swap = df.feature1 > df.feature2
    df.loc[swap, ["feature1", "feature2"]] = \
        df.loc[swap, ["feature2", "feature1"]].values
    return df.drop_duplicates(subset=["feature1", "feature2"]).reset_index(drop=True)


def load_genotypes(genotype_file, withdrawal_path, snps, prs_scripts_dir,
                   cache_file=None):
    """Pull just the needed SNP columns out of the .raw file, reusing
    prsInteractive's chunked reader so behaviour matches the rest of the
    pipeline (withdrawal removal, dtype handling, IID index).

    The training .raw is ~275 GB and a full scan takes ~10 minutes, so the
    extracted submatrix is cached when --cache_file is given."""
    if cache_file and os.path.exists(cache_file):
        print(f"loading cached genotypes from {cache_file} ...", flush=True)
        df = pd.read_parquet(cache_file)
        have = set(df.columns)
        if set(snps) - have - {"PHENOTYPE"}:
            print("cache does not cover every requested SNP; re-reading the .raw")
        else:
            return df

    sys.path.insert(0, prs_scripts_dir)
    from helper.download import get_dataset

    print(f"loading {len(snps)} SNP columns from {genotype_file} ...", flush=True)
    df = get_dataset(genotype_file, withdrawal_path, list(snps), use_chunking=True)
    if "PHENOTYPE" not in df.columns:
        raise SystemExit("genotype frame has no PHENOTYPE column")

    if cache_file:
        os.makedirs(os.path.dirname(cache_file) or ".", exist_ok=True)
        print(f"caching extracted genotypes to {cache_file} ...", flush=True)
        df.to_parquet(cache_file)
    return df


def build_design(geno, covar_file):
    """Merge covariates onto the genotype frame (indexed by IID) and return
    (X_covariates, y, genotype_frame) aligned on the same individuals."""
    covar = pd.read_csv(covar_file)
    # some covariate extracts spell it "Sex" rather than "SEX"
    if "SEX" not in covar.columns and "Sex" in covar.columns:
        covar = covar.rename(columns={"Sex": "SEX"})
    if covar["SEX"].dtype == object:
        covar["SEX"] = covar["SEX"].map(
            {"Male": 1, "Female": 0, "male": 1, "female": 0}).fillna(covar["SEX"])
    missing = [c for c in COVARIATES + ["IID"] if c not in covar.columns]
    if missing:
        raise SystemExit(f"covariate file missing columns: {missing}")
    covar = covar[["IID"] + COVARIATES].drop_duplicates(subset="IID")
    covar["IID"] = covar["IID"].astype(int)
    covar = covar.set_index("IID")

    y = geno["PHENOTYPE"].astype(float)
    if np.nanmax(y.values) == 2:            # plink 1/2 case-control coding
        y = y - 1
    geno = geno.drop(columns=["PHENOTYPE"])

    common = geno.index.intersection(covar.index)
    print(f"individuals: genotypes {len(geno)}, covariates {len(covar)}, "
          f"overlap {len(common)}")
    geno = geno.loc[common]
    y = y.loc[common]
    C = covar.loc[common, COVARIATES].astype(float)

    ok = C.notna().all(axis=1) & y.notna()
    dropped = int((~ok).sum())
    if dropped:
        print(f"dropping {dropped} individuals with missing covariates/phenotype")
    C, y, geno = C[ok], y[ok], geno[ok]

    # Mean-impute dosages once, up front. Doing it here rather than per pair
    # keeps the worker functions read-only, which matters because joblib
    # memory-maps large arrays to share them and those views cannot be written.
    geno = geno.astype(np.float32)
    n_missing = int(geno.isna().to_numpy().sum())
    if n_missing:
        print(f"mean-imputing {n_missing} missing dosages")
        geno = geno.fillna(geno.mean())
    return C, y, geno


# ------------------------------------------------------------------- the model

def fit_pair(g1, g2, C, y, min_double_carriers):
    """Fit logit(y) ~ C + b1*g1 + b2*g2 + J*(g1*g2).

    Returns a dict with the three effect estimates and diagnostics. J is the
    departure from additivity; b1/b2 are the constituent marginal effects that
    the standalone GxG models never estimated."""
    import statsmodels.api as sm

    n_double = int(((g1 > 0) & (g2 > 0)).sum())
    # Centre the dosages before forming the product. J, its SE and p-value are
    # invariant to this (the model is only reparameterised), but it takes the
    # collinearity between the product and its constituents from a median |r|
    # of 0.63 to 0.004, so b1 and b2 become the effect of each locus at the
    # mean of the other -- i.e. marginal effects -- rather than the effect when
    # the partner locus is zero, which is what the raw parameterisation gives.
    g1 = g1 - g1.mean()
    g2 = g2 - g2.mean()
    prod = g1 * g2
    out = {
        "n": int(len(y)),
        "n_double_carriers": n_double,
        "maf1": float(np.nanmean(g1) / 2),
        "maf2": float(np.nanmean(g2) / 2),
        "converged": False,
        "note": "",
    }
    for k in ("beta1", "se1", "p1", "beta2", "se2", "p2", "J", "se_J", "p_J"):
        out[k] = np.nan

    if n_double < min_double_carriers:
        out["note"] = f"below min_double_carriers ({n_double})"
        return out
    if prod.std() == 0 or g1.std() == 0 or g2.std() == 0:
        out["note"] = "no variance in a term"
        return out

    X = np.column_stack([np.ones(len(y)), C.values, g1, g2, prod])
    try:
        res = sm.Logit(y.values, X).fit(disp=0, maxiter=100)
    except Exception as exc:                       # separation, singular matrix
        out["note"] = f"fit failed: {type(exc).__name__}"
        return out

    # ones + covariates occupy the leading columns; g1, g2, product are last 3
    i1, i2, ij = X.shape[1] - 3, X.shape[1] - 2, X.shape[1] - 1
    out.update({
        "beta1": res.params[i1], "se1": res.bse[i1], "p1": res.pvalues[i1],
        "beta2": res.params[i2], "se2": res.bse[i2], "p2": res.pvalues[i2],
        "J": res.params[ij], "se_J": res.bse[ij], "p_J": res.pvalues[ij],
        "converged": bool(res.mle_retvals.get("converged", False)),
    })
    return out


def benjamini_hochberg(p):
    """BH q-values; NaN p-values pass through as NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty(m)
    tmp[order] = np.clip(adj, 0, 1)
    q[ok] = tmp
    return q


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs_file", required=True)
    ap.add_argument("--genotype_file", required=True,
                    help="training .raw file (the model is fit on training data)")
    ap.add_argument("--covar_file", required=True,
                    help="CSV with IID, age, SEX, PC1..PC10")
    ap.add_argument("--withdrawal_path", required=True)
    ap.add_argument("--out_dir", required=True,
                    help="new directory for results; existing files are not touched")
    ap.add_argument("--prs_scripts_dir",
                    default="/Users/kerimulterer/prsInteractive/scripts",
                    help="for the helper.download chunked .raw reader")
    ap.add_argument("--min_double_carriers", type=int, default=100,
                    help="pairs with fewer double carriers are not fit (default 100)")
    ap.add_argument("--q_threshold", type=float, default=0.05,
                    help="BH q-value cutoff on J (default 0.05)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--cache_file", default=None,
                    help="parquet path to cache the extracted genotype submatrix; "
                         "scanning the 275 GB training .raw takes ~10 minutes, so "
                         "this makes re-runs and threshold sweeps cheap")
    ap.add_argument("--limit", type=int, default=None,
                    help="test only the first N pairs (for a quick check)")
    args = ap.parse_args()

    real_out = os.path.realpath(args.out_dir)
    if os.path.realpath("/Users/kerimulterer/prsInteractive/results") in real_out:
        sys.exit("refusing to write inside prsInteractive/results; "
                 "choose an --out_dir outside it")
    os.makedirs(args.out_dir, exist_ok=True)

    pairs = load_pairs(args.pairs_file)
    if args.limit:
        pairs = pairs.head(args.limit)
    snps = sorted(set(pairs.feature1) | set(pairs.feature2))
    print(f"{len(pairs)} unique pairs over {len(snps)} distinct SNPs")

    geno = load_genotypes(args.genotype_file, args.withdrawal_path, snps,
                          args.prs_scripts_dir, args.cache_file)
    available = set(geno.columns)
    missing = [s for s in snps if s not in available]
    if missing:
        print(f"WARNING: {len(missing)} SNPs absent from the genotype file "
              f"(e.g. {missing[:5]}); their pairs are skipped")

    C, y, geno = build_design(geno, args.covar_file)
    print(f"fitting on n={len(y)} individuals, {int(y.sum())} cases "
          f"({y.mean():.2%})", flush=True)

    testable = pairs[pairs.feature1.isin(available) & pairs.feature2.isin(available)]
    skipped = pairs[~pairs.index.isin(testable.index)]
    print(f"testable pairs: {len(testable)} ({len(skipped)} skipped for missing SNPs)")

    from joblib import Parallel, delayed
    st = time.time()

    # dosages are already imputed and never written to, so workers share them
    def work(f1, f2):
        r = fit_pair(np.asarray(geno[f1], dtype=float),
                     np.asarray(geno[f2], dtype=float),
                     C, y, args.min_double_carriers)
        r["feature1"], r["feature2"] = f1, f2
        return r

    rows = Parallel(n_jobs=args.jobs, verbose=5)(
        delayed(work)(r.feature1, r.feature2) for r in testable.itertuples(index=False))
    print(f"fitting took {(time.time()-st)/60:.1f} minutes")

    res = pd.DataFrame(rows)
    res["q_J"] = benjamini_hochberg(res["p_J"].values)

    # carry the original pipeline coefficient through for comparison
    if "mean_coef" in pairs.columns:
        res = res.merge(pairs[["feature1", "feature2", "mean_coef"]],
                        on=["feature1", "feature2"], how="left")
        res = res.rename(columns={"mean_coef": "pipeline_pair_coef"})

    front = ["feature1", "feature2", "n", "n_double_carriers", "maf1", "maf2",
             "beta1", "se1", "p1", "beta2", "se2", "p2",
             "J", "se_J", "p_J", "q_J", "converged", "note"]
    res = res[[c for c in front if c in res.columns] +
              [c for c in res.columns if c not in front]]

    keep = (res.converged & res.q_J.notna() & (res.q_J < args.q_threshold))
    res["non_additive"] = keep

    res_path = os.path.join(args.out_dir, "gxg_joint_model_results.csv")
    res.to_csv(res_path, index=False)
    res[keep].to_csv(os.path.join(args.out_dir, "gxg_pairs.non_additive.csv"),
                     index=False)

    excl = res[~keep].copy()
    excl["reason"] = np.where(
        excl.note != "", excl.note,
        np.where(~excl.converged, "model did not converge",
                 f"J not significant at q<{args.q_threshold}"))
    excl.to_csv(os.path.join(args.out_dir, "gxg_pairs.excluded.csv"), index=False)

    n_lowcount = int((res.note.str.startswith("below min_double_carriers")).sum())
    n_failed = int(((res.note.str.startswith("fit failed")) |
                    (~res.converged & (res.note == ""))).sum())
    summary = (
        f"filter_non_additive_gxg_features.py\n"
        f"pairs file              : {args.pairs_file}\n"
        f"genotypes               : {args.genotype_file}\n"
        f"individuals fit         : {len(y)} ({int(y.sum())} cases)\n"
        f"covariates              : {', '.join(COVARIATES)}\n"
        f"min double carriers     : {args.min_double_carriers}\n"
        f"BH q threshold on J     : {args.q_threshold}\n"
        f"\n"
        f"pairs in input          : {len(pairs)}\n"
        f"skipped, SNP missing    : {len(skipped)}\n"
        f"skipped, low carriers   : {n_lowcount}\n"
        f"model failures          : {n_failed}\n"
        f"tested with a valid J   : {int(res.q_J.notna().sum())}\n"
        f"NON-ADDITIVE (retained) : {int(keep.sum())}\n"
        f"dropped, J not signif.  : {int((res.q_J.notna() & ~keep).sum())}\n"
    )
    with open(os.path.join(args.out_dir, "run_summary.txt"), "w") as f:
        f.write(summary)
    print("\n" + summary)
    print(f"results -> {res_path}")


if __name__ == "__main__":
    main()
