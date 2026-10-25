#!/usr/bin/env python3
"""
holdout_auc_table.py

Computes AUC with DeLong 95% confidence intervals on the holdout split for the
retained ePRS models, alongside the published validation values, for Table 2.

The holdout cohort was excluded from all model development (Methods 4.2.2),
whereas the validation split was used for SHAP-based feature selection
(Methods 4.5.4) and the model redundancy decisions (Methods 4.7). Reporting
both closes the question of how much optimism that selection introduced.

DeLong's method is used for the interval so it is consistent with the DeLong
test already used for the pairwise AUC comparisons in Table 1.

Outputs a plain table, a CSV, and LaTeX rows ready to paste.
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

# prsInteractive model name -> manuscript label
RETAINED = {
    "type2Diabetes": [("main", "PRS_G"), ("epi_summed", "PRS_{G+G}"),
                      ("epi_product", "PRS_{GxG}"), ("cardio_summed", "PRS_{G+GxE}"),
                      ("cardio_product", "PRS_{GxGxE}")],
    "celiacDisease": [("main", "PRS_G"), ("epi_product", "PRS_{GxG}"),
                      ("cardio_summed", "PRS_{G+GxE}")],
}


def midrank(x):
    """Midranks, as required by DeLong's structural components."""
    order = np.argsort(x)
    s = x[order]
    n = len(x)
    r = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n - 1 and s[j + 1] == s[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    out = np.empty(n, dtype=float)
    out[order] = r
    return out


def delong_auc_ci(y, score, alpha=0.05):
    """AUC and DeLong confidence interval for a single classifier."""
    pos = score[y == 1]
    neg = score[y == 0]
    m, n = len(pos), len(neg)

    tx, ty, tz = midrank(pos), midrank(neg), midrank(np.concatenate([pos, neg]))
    auc = (tz[:m].sum() - m * (m + 1) / 2) / (m * n)
    v01 = (tz[:m] - tx) / n                 # structural components, cases
    v10 = 1.0 - (tz[m:] - ty) / m           # structural components, controls
    var = v01.var(ddof=1) / m + v10.var(ddof=1) / n
    se = np.sqrt(var)
    z = stats.norm.ppf(1 - alpha / 2)
    return auc, max(0.0, auc - z * se), min(1.0, auc + z * se), se


def score_path(scores_dir, split):
    return os.path.join(scores_dir, "combinedPRSGroups.csv" if split == "validation"
                        else f"combinedPRSGroups.{split}.csv")


def load(scores_dir, model, split):
    """Read one model's score from combinedPRSGroups, which is what the
    published AUCs were computed from (prsAUCDelongStats.R fits
    PHENOTYPE ~ prs and takes the AUC of the fitted probabilities; with a
    single predictor that is a monotonic transform, so it equals the AUC of
    the score itself).

    The per-model prsScores/*.prs.csv files are not always in step with this
    file. For celiacDisease cardio_summed they differ by 0.026 in AUC (0.656
    against 0.630), and only the combinedPRSGroups value reproduces Table 2.
    """
    d = pd.read_csv(score_path(scores_dir, split),
                    usecols=["IID", "PHENOTYPE", f"prs_{model}"])
    return d.rename(columns={f"prs_{model}": "prs"})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results_root", default="/Users/kerimulterer/prsInteractive/results")
    ap.add_argument("--phenotypes", default="type2Diabetes,celiacDisease")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rows = []

    for pheno in [p.strip() for p in args.phenotypes.split(",")]:
        scores_dir = os.path.join(args.results_root, pheno,
                                  "combinedAnalysis", "scores")
        if not os.path.isdir(scores_dir):
            print(f"skipping {pheno}: no scores directory")
            continue
        for model, label in RETAINED.get(pheno, []):
            rec = {"phenotype": pheno, "model": model, "label": label}
            for split in ("validation", "holdout"):
                path = score_path(scores_dir, split)
                if not os.path.exists(path):
                    rec[f"{split}_auc"] = np.nan
                    continue
                d = load(scores_dir, model, split)
                y = d["PHENOTYPE"].astype(float).to_numpy()
                if np.nanmax(y) == 2:
                    y = y - 1
                s = d["prs"].to_numpy(dtype=float)
                ok = np.isfinite(y) & np.isfinite(s)
                auc, lo, hi, se = delong_auc_ci(y[ok], s[ok])
                rec[f"{split}_n"] = int(ok.sum())
                rec[f"{split}_cases"] = int(y[ok].sum())
                rec[f"{split}_auc"] = auc
                rec[f"{split}_lo"] = lo
                rec[f"{split}_hi"] = hi
            rec["delta"] = rec.get("holdout_auc", np.nan) - rec.get("validation_auc", np.nan)
            rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.out_dir, "table2_holdout_auc.csv"), index=False)

    fmt = lambda a, l, h: f"{a:.2f} ({l:.2f}–{h:.2f})"
    for pheno, grp in out.groupby("phenotype", sort=False):
        v_n, v_c = int(grp.validation_n.iloc[0]), int(grp.validation_cases.iloc[0])
        h_n, h_c = int(grp.holdout_n.iloc[0]), int(grp.holdout_cases.iloc[0])
        print(f"\n{pheno}  validation n={v_n} ({v_c} cases) | "
              f"holdout n={h_n} ({h_c} cases)")
        print(f'{"model":14} {"validation AUC":>22} {"holdout AUC":>22} {"delta":>7}')
        for r in grp.itertuples():
            print(f"{r.label:14} {fmt(r.validation_auc, r.validation_lo, r.validation_hi):>22} "
                  f"{fmt(r.holdout_auc, r.holdout_lo, r.holdout_hi):>22} {r.delta:>+7.3f}")

    tex = [r"% Holdout AUC row for Table 2 -- paste under the existing AUC row.",
           r"% Holdout n and cases are given in the caption, not per column."]
    for pheno, grp in out.groupby("phenotype", sort=False):
        cells = " & ".join(fmt(r.validation_auc, r.validation_lo, r.validation_hi)
                           for r in grp.itertuples())
        hcells = " & ".join(fmt(r.holdout_auc, r.holdout_lo, r.holdout_hi)
                            for r in grp.itertuples())
        tex.append(f"% --- {pheno} ({', '.join(g.label for g in grp.itertuples())})")
        tex.append(rf"AUC, validation (95\% CI) & {cells} \\")
        tex.append(rf"AUC, holdout (95\% CI) & {hcells} \\")
    path = os.path.join(args.out_dir, "table2_holdout_auc.tex")
    with open(path, "w") as f:
        f.write("\n".join(tex) + "\n")
    print(f"\nLaTeX rows -> {path}")


if __name__ == "__main__":
    main()
