#!/usr/bin/env python3
"""
filter_non_additive_gxe_features.py

Departure-from-additivity test for the GxE (and GxGxE) features of the ePRS
pipeline — the environmental counterpart of filter_non_additive_gxg_features.py.

WHAT prsInteractive BUILDS
--------------------------
Per Methods 4.4.3, an environmental interaction feature is a single column
formed by multiplying a centred environmental measure by an encoded genetic
feature:

    G+GxE  :  Z = E * (g_j + g_l)        environment x additive burden
    GxGxE  :  Z = E * g_j * g_l          environment x epistatic product

and a two-locus pair may instead be a single SNP or an imputed HLA allele, in
which case the feature is simply E * g.

WHY A JOINT MODEL IS NEEDED
---------------------------
An interaction coefficient only means "departure from additivity" if every
lower-order term is in the same model. filter_non_additive_gen_env_features.py
compares a GxE coefficient against two things: the environmental main effect E,
and the coefficient of the genetic composite ("g_j,g_l"). For the product
encoding that leaves four lower-order terms unaccounted for — g_j, g_l, E*g_j
and E*g_l — so the three-way coefficient absorbs them. It is also a comparison
of magnitudes across separately fitted models, with no standard error and so no
test.

This script fits the full hierarchy in one model per feature:

  two-way   E * g          y ~ C + E + g + J*(E*g)
  summed    E * (g_j+g_l)  y ~ C + E + g_j + g_l + J*(E*(g_j+g_l))
  product   E * g_j * g_l  y ~ C + E + g_j + g_l + g_j*g_l
                                 + E*g_j + E*g_l + J*(E*g_j*g_l)

J is then the departure from additivity at the highest order, with a Wald
p-value and a BH q-value. The summed form is a *constrained* two-way
interaction (it forces E to interact equally with both loci), which is what the
pipeline's feature construction assumes; the product form is a genuine
three-way interaction and needs all six lower-order terms.

Every constituent is mean-centred before products are formed. As established
for GxG, centring leaves the highest-order coefficient and its p-value
unchanged while removing the collinearity that otherwise makes the lower-order
coefficients uninterpretable.

OUTPUTS (under --out_dir; nothing in prsInteractive is modified)
    gxe_joint_model_results.csv   every term tested, with all model terms,
                                  J, SE, p, q, carrier counts, convergence
    gxe_terms.non_additive.csv    terms passing
    gxe_terms.excluded.csv        terms failing, with a reason
    run_summary.txt
"""

import argparse
import glob
import os
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter(action="ignore")

COVARIATES = ["age", "SEX"] + [f"PC{i}" for i in range(1, 11)]
SNP_RE = re.compile(r"^(rs\d+|Affx-\S+)$")


def is_snp(name):
    return bool(SNP_RE.match(str(name)))


# ---------------------------------------------------------------- term loading

def load_terms(decoupled_dir, encoding, gxe_file=None):
    """Collect GxE terms for one encoding. Encoding is taken from the model
    directory name (*_product / *_summed), because it determines how the
    genetic composite was built and therefore which model is correct."""
    if gxe_file:
        df = pd.read_csv(gxe_file)
        df["model"] = os.path.basename(os.path.dirname(gxe_file))
        frames = [df]
    else:
        frames = []
        for path in sorted(glob.glob(os.path.join(decoupled_dir, "*", "GxE_features.csv"))):
            model = os.path.basename(os.path.dirname(path))
            if not model.endswith(f"_{encoding}"):
                continue
            d = pd.read_csv(path)
            d["model"] = model
            frames.append(d)
    if not frames:
        raise SystemExit(f"no GxE feature files found for encoding '{encoding}'")

    df = pd.concat(frames, ignore_index=True)
    df["feature2"] = df["feature2"].fillna("").astype(str)
    df["feature1"] = df["feature1"].astype(str)
    df["env"] = df["env"].astype(str)
    agg = (df.groupby(["env", "feature1", "feature2"], as_index=False)
             .agg(pipeline_coef=("coef", "mean"), n_models=("coef", "size")))
    return agg


# ---------------------------------------------------------------- data loading

def load_snp_dosages(genotype_file, withdrawal_path, snps, prs_scripts_dir, cache_file):
    if cache_file and os.path.exists(cache_file):
        print(f"loading cached genotypes from {cache_file} ...", flush=True)
        cached = pd.read_parquet(cache_file)
        missing = [s for s in snps if s not in cached.columns]
        if not missing:
            return cached
        print(f"cache is missing {len(missing)} SNPs; reading the .raw for the full set")

    sys.path.insert(0, prs_scripts_dir)
    from helper.download import get_dataset
    print(f"loading {len(snps)} SNP columns from {genotype_file} ...", flush=True)
    df = get_dataset(genotype_file, withdrawal_path, list(snps), use_chunking=True)
    if cache_file:
        os.makedirs(os.path.dirname(cache_file) or ".", exist_ok=True)
        print(f"caching to {cache_file} ...", flush=True)
        df.to_parquet(cache_file)
    return df


def load_hla(hla_file, alleles):
    if not alleles:
        return None
    print(f"loading {len(alleles)} HLA allelotypes from {hla_file} ...", flush=True)
    hla = pd.read_csv(hla_file, usecols=["Participant ID"] + list(alleles))
    hla = hla.rename(columns={"Participant ID": "IID"}).drop_duplicates("IID")
    hla["IID"] = hla["IID"].astype(int)
    return hla.set_index("IID")


# ------------------------------------------------------------------- the model

def centre(x):
    return x - np.nanmean(x)


def fit_term(E, g1, g2, C, y, encoding, min_carriers):
    """Fit the full hierarchy for one GxE term and return the highest-order
    interaction coefficient J plus every lower-order term."""
    import statsmodels.api as sm

    out = {"n": int(len(y)), "converged": False, "note": "", "n_carriers": 0}
    for k in ("beta_E", "se_E", "p_E", "beta_g1", "beta_g2", "beta_g1g2",
              "beta_Eg1", "beta_Eg2", "J", "se_J", "p_J"):
        out[k] = np.nan

    Ec = centre(E)
    g1c = centre(g1)

    if g2 is None:                                   # two-way: E * g
        out["n_carriers"] = int((g1 > 0).sum())
        cols = [Ec, g1c]
        names = ["beta_E", "beta_g1"]
        inter = Ec * g1c
    elif encoding == "summed":                       # E * (g_j + g_l)
        g2c = centre(g2)
        out["n_carriers"] = int(((g1 + g2) > 0).sum())
        cols = [Ec, g1c, g2c]
        names = ["beta_E", "beta_g1", "beta_g2"]
        inter = Ec * (g1c + g2c)
    else:                                            # E * g_j * g_l
        g2c = centre(g2)
        out["n_carriers"] = int(((g1 > 0) & (g2 > 0)).sum())
        cols = [Ec, g1c, g2c, g1c * g2c, Ec * g1c, Ec * g2c]
        names = ["beta_E", "beta_g1", "beta_g2", "beta_g1g2", "beta_Eg1", "beta_Eg2"]
        inter = Ec * g1c * g2c

    if out["n_carriers"] < min_carriers:
        out["note"] = f"below min_carriers ({out['n_carriers']})"
        return out
    if np.std(inter) == 0 or any(np.std(c) == 0 for c in cols):
        out["note"] = "no variance in a term"
        return out

    X = np.column_stack([np.ones(len(y)), C.values] + cols + [inter])
    try:
        res = sm.Logit(y.values, X).fit(disp=0, maxiter=200)
    except Exception as exc:
        out["note"] = f"fit failed: {type(exc).__name__}"
        return out

    base = 1 + C.shape[1]
    for i, nm in enumerate(names):
        out[nm] = res.params[base + i]
    out["J"] = res.params[-1]
    out["se_J"] = res.bse[-1]
    out["p_J"] = res.pvalues[-1]
    out["se_E"] = res.bse[base]
    out["p_E"] = res.pvalues[base]
    out["converged"] = bool(res.mle_retvals.get("converged", False))
    out["n_terms"] = len(names) + 1
    return out


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    m = len(ranked)
    adj = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    tmp = np.empty(m)
    tmp[order] = np.clip(adj, 0, 1)
    q[ok] = tmp
    return q


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--decoupled_dir", required=True)
    ap.add_argument("--gxe_file", default=None,
                    help="a single GxE_features.csv; otherwise all models matching "
                         "--encoding are pooled")
    ap.add_argument("--encoding", required=True, choices=["product", "summed"])
    ap.add_argument("--env_file", required=True,
                    help="allEnvironmentTraining.csv (standardized E values)")
    ap.add_argument("--genotype_file", required=True)
    ap.add_argument("--covar_file", required=True)
    ap.add_argument("--withdrawal_path", required=True)
    ap.add_argument("--hla_file", default=None,
                    help="participant_hla.csv, needed for terms whose genetic "
                         "constituent is an imputed HLA allelotype")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--prs_scripts_dir",
                    default="/Users/kerimulterer/prsInteractive/scripts")
    ap.add_argument("--cache_file", default=None)
    ap.add_argument("--min_carriers", type=int, default=2000)
    ap.add_argument("--q_threshold", type=float, default=0.05)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    if os.path.realpath("/Users/kerimulterer/prsInteractive/results") in \
            os.path.realpath(args.out_dir):
        sys.exit("refusing to write inside prsInteractive/results")
    os.makedirs(args.out_dir, exist_ok=True)

    terms = load_terms(args.decoupled_dir, args.encoding, args.gxe_file)
    if args.limit:
        terms = terms.head(args.limit)
    print(f"{len(terms)} unique {args.encoding} GxE terms")

    gen = set(terms.feature1) | {f for f in terms.feature2 if f}
    snps = sorted(g for g in gen if is_snp(g))
    hlas = sorted(g for g in gen if g and not is_snp(g))
    print(f"genetic constituents: {len(snps)} SNPs, {len(hlas)} HLA allelotypes")

    geno = load_snp_dosages(args.genotype_file, args.withdrawal_path, snps,
                            args.prs_scripts_dir, args.cache_file)
    y = geno["PHENOTYPE"].astype(float)
    if np.nanmax(y.values) == 2:
        y = y - 1
    geno = geno.drop(columns=["PHENOTYPE"])

    hla = load_hla(args.hla_file, hlas) if (hlas and args.hla_file) else None
    if hlas and hla is None:
        print(f"WARNING: {len(hlas)} HLA constituents but no --hla_file; "
              f"their terms will be skipped")

    env = pd.read_csv(args.env_file)
    env["IID"] = env["IID"].astype(int)
    env = env.drop_duplicates("IID").set_index("IID")

    covar = pd.read_csv(args.covar_file)
    if "SEX" not in covar.columns and "Sex" in covar.columns:
        covar = covar.rename(columns={"Sex": "SEX"})
    covar = covar[["IID"] + COVARIATES].drop_duplicates("IID")
    covar["IID"] = covar["IID"].astype(int)
    covar = covar.set_index("IID")

    idx = geno.index.intersection(covar.index).intersection(env.index)
    if hla is not None:
        idx = idx.intersection(hla.index)
    print(f"individuals with genotypes, covariates, environment"
          f"{' and HLA' if hla is not None else ''}: {len(idx)}")

    C = covar.loc[idx, COVARIATES].astype(float)
    ok = C.notna().all(axis=1) & y.loc[idx].notna()
    idx = idx[ok.values]
    C = C.loc[idx]
    y = y.loc[idx]
    geno = geno.loc[idx].astype(np.float32)
    geno = geno.fillna(geno.mean())
    env = env.loc[idx]
    # allEnvironment*.csv still carries missing values (e.g. ~13% for Albumin),
    # and statsmodels rejects NaN outright. Mean-impute to match how the
    # pipeline's SimpleImputer treats the same columns downstream.
    n_env_missing = int(env.isna().to_numpy().sum())
    if n_env_missing:
        print(f"mean-imputing {n_env_missing} missing environmental values")
        env = env.fillna(env.mean())
    if hla is not None:
        hla = hla.loc[idx].astype(np.float32)
        hla = hla.fillna(hla.mean())
    print(f"fitting on n={len(y)}, {int(y.sum())} cases ({y.mean():.2%})", flush=True)

    env_avail = set(env.columns)
    def genetic_ok(name):
        if not name:
            return True
        return name in geno.columns or (hla is not None and name in hla.columns)

    testable = terms[terms.env.isin(env_avail)
                     & terms.feature1.map(genetic_ok)
                     & terms.feature2.map(genetic_ok)]
    skipped = len(terms) - len(testable)
    if skipped:
        print(f"{skipped} terms skipped: environmental or genetic constituent unavailable")
    print(f"testable terms: {len(testable)}", flush=True)

    def col(name):
        if name in geno.columns:
            return np.asarray(geno[name], dtype=float)
        return np.asarray(hla[name], dtype=float)

    from joblib import Parallel, delayed

    def work(env_name, f1, f2):
        r = fit_term(np.asarray(env[env_name], dtype=float),
                     col(f1), col(f2) if f2 else None,
                     C, y, args.encoding, args.min_carriers)
        r["env"], r["feature1"], r["feature2"] = env_name, f1, f2
        r["order"] = "two-way" if not f2 else ("three-way" if args.encoding == "product"
                                               else "summed-pair")
        return r

    st = time.time()
    rows = Parallel(n_jobs=args.jobs, verbose=5)(
        delayed(work)(t.env, t.feature1, t.feature2)
        for t in testable.itertuples(index=False))
    print(f"fitting took {(time.time()-st)/60:.1f} minutes")

    res = pd.DataFrame(rows)
    res["q_J"] = benjamini_hochberg(res["p_J"].values)
    res = res.merge(testable[["env", "feature1", "feature2", "pipeline_coef"]],
                    on=["env", "feature1", "feature2"], how="left")
    keep = res.converged & res.q_J.notna() & (res.q_J < args.q_threshold)
    res["non_additive"] = keep

    front = ["env", "feature1", "feature2", "order", "n", "n_carriers",
             "beta_E", "se_E", "p_E", "beta_g1", "beta_g2", "beta_g1g2",
             "beta_Eg1", "beta_Eg2", "J", "se_J", "p_J", "q_J",
             "pipeline_coef", "converged", "note", "non_additive"]
    res = res[[c for c in front if c in res.columns] +
              [c for c in res.columns if c not in front]]

    res.to_csv(os.path.join(args.out_dir, "gxe_joint_model_results.csv"), index=False)
    res[keep].to_csv(os.path.join(args.out_dir, "gxe_terms.non_additive.csv"), index=False)
    ex = res[~keep].copy()
    ex["reason"] = np.where(ex.note != "", ex.note,
                            np.where(~ex.converged, "did not converge",
                                     f"J not significant at q<{args.q_threshold}"))
    ex.to_csv(os.path.join(args.out_dir, "gxe_terms.excluded.csv"), index=False)

    n_low = int(res.note.str.startswith("below min_carriers").sum())
    summary = (
        f"filter_non_additive_gxe_features.py\n"
        f"encoding                : {args.encoding}\n"
        f"decoupled dir           : {args.decoupled_dir}\n"
        f"individuals fit         : {len(y)} ({int(y.sum())} cases)\n"
        f"covariates              : {', '.join(COVARIATES)}\n"
        f"min carriers            : {args.min_carriers}\n"
        f"BH q threshold on J     : {args.q_threshold}\n\n"
        f"terms in input          : {len(terms)}\n"
        f"skipped, constituent NA : {skipped}\n"
        f"skipped, low carriers   : {n_low}\n"
        f"tested with a valid J   : {int(res.q_J.notna().sum())}\n"
        f"NON-ADDITIVE (retained) : {int(keep.sum())}\n"
        f"dropped, J not signif.  : {int((res.q_J.notna() & ~keep).sum())}\n")
    with open(os.path.join(args.out_dir, "run_summary.txt"), "w") as f:
        f.write(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
