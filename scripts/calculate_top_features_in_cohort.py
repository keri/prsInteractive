#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
import re

from sklearn.metrics import roc_auc_score, make_scorer, roc_curve, auc
from sklearn.impute import SimpleImputer
import warnings

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from scipy import stats as scipy_stats
import os

from helper.download import *
from helper.data_wrangling import *
from helper.draw_plots import *
from helper.calculate_shap_values import *

warnings.filterwarnings("ignore", message="Bad circle positioning")


# =============================================================================
#  FEATURE FILTERING STRATEGIES
#  Switch between these by passing --filter_strategy to main()
#
#  1. "strict"       – original: positive in matching cohort AND ALL others negative
#  2. "majority"     – positive in matching cohort AND negative in >=majority_threshold fraction of others
#  3. "differential" – matching cohort z-score is highest AND delta vs mean-of-others > delta_threshold
#  4. "tiered"       – runs all three and labels each feature with the strictest tier it qualifies for
# =============================================================================


def find_opposite_sign_features(df, feature_sets, threshold):
    """
    ORIGINAL (strict) strategy.
    Positive z_score in matching cohort AND ALL other cohorts negative.
    """
    rows = []
    for f in feature_sets:
        df_feature_set = df[df['features_used'] == f]
        df_matching = df_feature_set[(df_feature_set['cohort'] == f) & (df_feature_set['z_score'] > threshold)]
        df_other = df_feature_set[(df_feature_set['cohort'] != f) & (df_feature_set['feature'].isin(df_matching.feature.tolist()))]

        for feature in df_matching['feature'].unique():
            matching_row = df_matching[df_matching['feature'] == feature]
            if len(matching_row) == 0:
                continue

            matching_z = matching_row['z_score'].values[0]
            matching_coef = matching_row['coefs'].values[0]
            other_z = df_other[df_other['feature'] == feature]['z_score'].values

            if len(other_z) == 0:
                continue

            if matching_z > 0 and all(z < 0 for z in other_z):
                rows.append({
                    'feature_set': f,
                    'feature': feature,
                    'matching_z_score': matching_z,
                    'n_other_cohorts': len(other_z),
                    'coefs': matching_coef,
                    'odds_ratio': np.exp(matching_coef)
                })

    return pd.DataFrame(rows)


def find_majority_sign_features(df, feature_sets, threshold, majority_threshold=0.6):
    """
    RELAXED majority strategy.
    Positive z_score in matching cohort AND negative in >= majority_threshold fraction of other cohorts.

    majority_threshold=0.6 means "negative in at least 60% of other cohorts".
    This recovers features that are specific to a cohort even if 1-2 other cohorts
    happen to be weakly positive.

    Parameters
    ----------
    majority_threshold : float, 0–1
        Fraction of other cohorts that must have negative z-score.
        0.5 = >half, 0.75 = three-quarters, 1.0 = all (same as strict).
    """
    rows = []
    for f in feature_sets:
        df_feature_set = df[df['features_used'] == f]
        df_matching = df_feature_set[(df_feature_set['cohort'] == f) & (df_feature_set['z_score'] > threshold)]
        df_other = df_feature_set[(df_feature_set['cohort'] != f) & (df_feature_set['feature'].isin(df_matching.feature.tolist()))]

        for feature in df_matching['feature'].unique():
            matching_row = df_matching[df_matching['feature'] == feature]
            if len(matching_row) == 0:
                continue

            matching_z = matching_row['z_score'].values[0]
            matching_coef = matching_row['coefs'].values[0]
            other_rows = df_other[df_other['feature'] == feature]
            other_z = other_rows['z_score'].values

            if len(other_z) == 0:
                continue

            frac_negative = np.sum(other_z < 0) / len(other_z)
            if frac_negative == np.nan:
              print('matching_row = ',matching_row)
              print('other z = ',other_z)

            if matching_z > 0 and frac_negative >= majority_threshold:
                rows.append({
                    'feature_set': f,
                    'feature': feature,
                    'matching_z_score': matching_z,
                    'n_other_cohorts': len(other_z),
                    'frac_other_negative': frac_negative,
                    'coefs': matching_coef,
                    'odds_ratio': np.exp(matching_coef)
                })

    return pd.DataFrame(rows)


def find_differential_features(df, feature_sets, threshold, delta_threshold=0.5):
    """
    DIFFERENTIAL enrichment strategy.
    Keeps features where:
      - z_score in matching cohort > threshold, AND
      - matching z_score is the HIGHEST across all cohorts for that feature set, AND
      - (matching_z - mean_of_others) >= delta_threshold

    This does NOT require other cohorts to be negative — just that the matching
    cohort stands out clearly above the rest.  Useful when all cohorts trend
    positive but one cohort is distinctly higher.

    Parameters
    ----------
    delta_threshold : float
        Minimum required difference between matching z and mean of other cohorts.
    """
    rows = []
    for f in feature_sets:
        df_feature_set = df[df['features_used'] == f]
        df_matching = df_feature_set[(df_feature_set['cohort'] == f) & (df_feature_set['z_score'] > threshold)]
        df_other = df_feature_set[(df_feature_set['cohort'] != f) & (df_feature_set['feature'].isin(df_matching.feature.tolist()))]

        for feature in df_matching['feature'].unique():
            matching_row = df_matching[df_matching['feature'] == feature]
            if len(matching_row) == 0:
                continue

            matching_z = matching_row['z_score'].values[0]
            matching_coef = matching_row['coefs'].values[0]
            other_z = df_other[df_other['feature'] == feature]['z_score'].values

            if len(other_z) == 0:
                continue

            mean_other_z = np.mean(other_z)
            delta = matching_z - mean_other_z

            # Must be highest z across all cohorts AND clearly elevated
            is_highest = matching_z > np.max(other_z)

            if is_highest and delta >= delta_threshold:
                rows.append({
                    'feature_set': f,
                    'feature': feature,
                    'matching_z_score': matching_z,
                    'mean_other_z': mean_other_z,
                    'delta_z': delta,
                    'n_other_cohorts': len(other_z),
                    'coefs': matching_coef,
                    'odds_ratio': np.exp(matching_coef)
                })

    return pd.DataFrame(rows)


def find_tiered_features(df, feature_sets, threshold,
                         majority_threshold=0.6, delta_threshold=0.5):
    """
    TIERED strategy — runs all three filters and assigns each surviving feature
    the strictest tier it qualifies for:

      Tier 1 (most specific)  : positive in matching, ALL others negative (original)
      Tier 2 (moderate)       : positive in matching, majority of others negative
      Tier 3 (least specific) : matching is highest AND delta > delta_threshold

    This lets you use tier as a confidence label rather than binary keep/drop.
    """
    strict_df = find_opposite_sign_features(df, feature_sets, threshold)
    majority_df = find_majority_sign_features(df, feature_sets, threshold, majority_threshold)
    differential_df = find_differential_features(df, feature_sets, threshold, delta_threshold)

    # Assign tiers — start permissive, then override with stricter
    tier_rows = []

    all_features = set()
    if not differential_df.empty:
        all_features |= set(zip(differential_df['feature_set'], differential_df['feature']))
    if not majority_df.empty:
        all_features |= set(zip(majority_df['feature_set'], majority_df['feature']))
    if not strict_df.empty:
        all_features |= set(zip(strict_df['feature_set'], strict_df['feature']))

    for (fset, feat) in all_features:
        # Find the source row — prefer strict > majority > differential
        source_df = None
        tier = None

        if not strict_df.empty and len(strict_df[(strict_df['feature_set'] == fset) & (strict_df['feature'] == feat)]) > 0:
            source_df = strict_df
            tier = 1
        elif not majority_df.empty and len(majority_df[(majority_df['feature_set'] == fset) & (majority_df['feature'] == feat)]) > 0:
            source_df = majority_df
            tier = 2
        elif not differential_df.empty and len(differential_df[(differential_df['feature_set'] == fset) & (differential_df['feature'] == feat)]) > 0:
            source_df = differential_df
            tier = 3

        if source_df is None:
            continue

        row = source_df[(source_df['feature_set'] == fset) & (source_df['feature'] == feat)].iloc[0].to_dict()
        row['specificity_tier'] = tier
        tier_rows.append(row)

    result = pd.DataFrame(tier_rows)
    if not result.empty:
        result = result.sort_values(['specificity_tier', 'odds_ratio'], ascending=[True, False])
    return result


# =============================================================================
#  Utility: per-cohort feature count diagnostic
# =============================================================================

def print_feature_diagnostics(df, feature_sets, threshold):
    """
    Before cross-cohort filtering, print how many features pass per cohort.
    Helps diagnose whether sparsity is upstream (per-cohort) or at intersection.
    """
    print("\n===== PER-COHORT FEATURE DIAGNOSTICS =====")
    for f in feature_sets:
        df_fs = df[df['features_used'] == f]
        total = len(df_fs['feature'].unique())
        n_important = len(df_fs[df_fs['feature_importance'] == 1]['feature'].unique())
        n_above_threshold = len(df_fs[(df_fs['cohort'] == f) & (df_fs['z_score'] > threshold)]['feature'].unique())
        print(f"  {f}: total={total}  flagged_important={n_important}  z>{threshold}_in_matching={n_above_threshold}")
    print("===========================================\n")


# =============================================================================
#  Unchanged helpers below
# =============================================================================

# =============================================================================
#  EXCLUSIVE HIGH-CASES: IDENTIFICATION AND BOOTSTRAP GBM-SHAP
# =============================================================================

def detect_bin_thresholds(groupedDf, bin_cols):
    """
    Infer high-risk / low-risk thresholds from the actual bin scale present in
    the data, preserving the same top-20% / bottom-20% split regardless of
    whether bins run 1–10 (deciles) or 1–1000 (quantiles).

      1–10  scale : high > 8   (= 0.8 × 10),  low < 3   (= 0.2 × 10 + 1)
      1–1000 scale: high > 800 (= 0.8 × 1000), low < 201 (= 0.2 × 1000 + 1)
    """
    max_bin = float(groupedDf[bin_cols].max().max())
    high_threshold = 0.8 * max_bin
    low_threshold  = 0.2 * max_bin + 1
    print(f'  Bin scale: 1–{max_bin:.0f}  →  '
          f'high risk bin > {high_threshold:.0f}  |  '
          f'low risk bin < {low_threshold:.0f}')
    return high_threshold, low_threshold


def identify_exclusive_high_cases(groupedDf, bin_cols):
    """
    For each cohort, identify individuals who are exclusively high-risk in that
    cohort but NOT high-risk in any other cohort, restricted to cases
    (PHENOTYPE == 2).  Controls are low-risk within the same cohort
    (PHENOTYPE == 1).

    Thresholds are derived proportionally from the bin scale (top 20% / bottom
    20%) so the function is correct for both the 1–10 training file and the
    1–1000 holdout file.

    Parameters
    ----------
    groupedDf : DataFrame indexed by IID with bin_* and PHENOTYPE columns
    bin_cols  : list of bin_* column names to consider

    Returns
    -------
    dict  {cohort: {'case_index': Index, 'control_index': Index}}
    """
    high_thresh, low_thresh = detect_bin_thresholds(groupedDf, bin_cols)
    result = {}
    cohort_names = [c.replace('bin_', '') for c in bin_cols]
    for col, cohort in zip(bin_cols, cohort_names):
        is_high_this  = groupedDf[col] > high_thresh
        other_cols    = [c for c in bin_cols if c != col]
        is_high_other = (
            groupedDf[other_cols].gt(high_thresh).any(axis=1)
            if other_cols
            else pd.Series(False, index=groupedDf.index)
        )
        case_idx    = groupedDf[
            is_high_this & ~is_high_other & (groupedDf['PHENOTYPE'] == 2)
        ].index
        control_idx = groupedDf[
            (groupedDf[col] < low_thresh) & (groupedDf['PHENOTYPE'] == 1)
        ].index
        result[cohort] = {'case_index': case_idx, 'control_index': control_idx}
        print(f'  ExclusiveHigh {cohort}: {len(case_idx):>5} exclusive cases  '
              f'| {len(control_idx):>6} within-cohort controls')
    return result


def run_exclusive_high_bootstrap(
        cohort, case_index, control_index,
        feature_df, featuresDf, feature_str,
        n_iterations, z_threshold,
        figPath, featureImportanceFile,
        rng_seed=42):
    """
    Train a Histogram Gradient Boosted classifier (exclusive cases vs balanced
    controls) n_iterations times, bootstrapping a fresh balanced control sample
    each iteration.  SHAP values are extracted from each model and the per-feature
    mean SHAP across all case samples is accumulated.  After all iterations the
    accumulated means are averaged and z-scored; features whose |z| exceeds
    z_threshold are flagged as important.

    Design rationale — why bootstrap balanced controls?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Exclusive cases are by definition a small, biologically coherent group.
    Controls are plentiful and heterogeneous.  Training on the full imbalanced
    dataset inflates the majority class and suppresses the minority signal.
    Bootstrapping n balanced subsets of controls and averaging SHAP z-scores
    across iterations gives a stable, low-variance importance estimate without
    any single control sample composition dominating the result.

    Why gradient boosted trees rather than random forest?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GBMs build trees sequentially to correct residuals, making them more
    sensitive to subtle, non-linear interactions in small exclusive cohorts.
    HistGradientBoostingClassifier is the sklearn native implementation —
    fast, supports native missing-value handling, and is fully SHAP-compatible
    via shap.TreeExplainer.

    Parameters
    ----------
    cohort              : str, cohort label (e.g. 'G', 'GxG')
    case_index          : pandas Index of exclusive high-risk case IIDs
    control_index       : pandas Index of within-cohort low-risk control IIDs
    feature_df          : DataFrame indexed by IID with feature columns
                          (all features for this model, IID-indexed after
                           assign_groups has been called)
    featuresDf          : DataFrame with 'feature' and 'features_used' columns
                          (already renamed from 'model')
    feature_str         : str, model name identifying which rows of featuresDf
                          to use and which output label to assign
    n_iterations        : int, number of bootstrap iterations
    z_threshold         : float, |z_score| threshold to flag feature_importance=1
    figPath             : str, directory for saving plots
    featureImportanceFile : str, path to the shared output CSV (appended to)
    rng_seed            : int, base RNG seed (each iteration offset by index)
    """
    # ---- Feature columns for this model -------------------------------------
    model_features = featuresDf[
        featuresDf['features_used'] == feature_str
    ]['feature'].tolist()
    feat_cols = [f for f in model_features if f in feature_df.columns]

    if not feat_cols:
        print(f'    [SKIP] No features for {cohort} × {feature_str}')
        return

    # ---- Intersect IIDs with feature_df index --------------------------------
    case_ids    = [i for i in case_index    if i in feature_df.index]
    control_ids = [i for i in control_index if i in feature_df.index]
    n_cases     = len(case_ids)

    if n_cases < 5:
        print(f'    [SKIP] {cohort} × {feature_str}: only {n_cases} exclusive cases')
        return
    if len(control_ids) < 5:
        print(f'    [SKIP] {cohort} × {feature_str}: only {len(control_ids)} controls')
        return

    X_cases         = feature_df.loc[case_ids,    feat_cols]
    X_controls_pool = feature_df.loc[control_ids, feat_cols]
    y_cases         = np.ones(n_cases,  dtype=int)

    replace = len(control_ids) < n_cases     # sample with replacement only if pool < cases
    if replace:
        print(f'    [WARN] {cohort}: control pool ({len(control_ids)}) < n_cases '
              f'({n_cases}); sampling with replacement')

    rng = np.random.default_rng(rng_seed)

    # ---- Bootstrap loop -----------------------------------------------------
    # Accumulate sum of per-feature mean-SHAP-over-cases for each iteration
    shap_sum      = np.zeros(len(feat_cols), dtype=float)
    n_successful  = 0

    for it in range(n_iterations):
        ctrl_pos    = rng.choice(len(control_ids), size=n_cases, replace=replace)
        X_ctrl_iter = X_controls_pool.iloc[ctrl_pos]
        y_ctrl_iter = np.zeros(n_cases, dtype=int)

        X_iter = pd.concat([X_cases, X_ctrl_iter])
        y_iter = np.concatenate([y_cases, y_ctrl_iter])

        try:
            gbm = HistGradientBoostingClassifier(
                max_iter=200,
                learning_rate=0.05,
                max_depth=4,
                random_state=int(rng.integers(0, 2**31)),
            )
            gbm.fit(X_iter, y_iter)

            # SHAP on cases only — we want the attribution that drives case
            # predictions, not the average across both classes.
            explainer  = shap.TreeExplainer(gbm)
            shap_vals  = explainer.shap_values(X_cases)

            # HistGBM binary: shap_values is a 2D array (n_samples × n_features)
            # for the positive class; no list wrapping needed.
            if isinstance(shap_vals, list):
                shap_vals = shap_vals[1]

            shap_sum     += np.mean(shap_vals, axis=0)
            n_successful += 1

        except Exception as exc:
            if it < 3:
                print(f'    [WARN] iteration {it} failed ({cohort} × {feature_str}): {exc}')
            continue

    if n_successful == 0:
        print(f'    [ERROR] All {n_iterations} iterations failed for '
              f'{cohort} × {feature_str}')
        return

    # ---- Average and z-score ------------------------------------------------
    avg_mean_shap = pd.Series(shap_sum / n_successful, index=feat_cols)

    # Drop features whose averaged SHAP is identically zero across all iterations
    avg_mean_shap = avg_mean_shap[avg_mean_shap != 0]
    if avg_mean_shap.empty:
        print(f'    [SKIP] {cohort} × {feature_str}: all averaged SHAP values zero')
        return

    z_scores = pd.Series(
        scipy_stats.zscore(avg_mean_shap.values),
        index=avg_mean_shap.index,
        name='z_score',
    )

    # ---- Build output row (same schema as the existing RF loop) -------------
    featureImportance = pd.DataFrame({
        'feature'   : avg_mean_shap.index,
        'mean_shap' : avg_mean_shap.values,
        'z_score'   : z_scores.values,
    })
    featureImportance['feature_importance'] = (
        featureImportance['z_score'].abs() > z_threshold
    ).astype(int)
    featureImportance['cohort']             = cohort
    featureImportance['features_used']      = feature_str
    featureImportance['training_data_used'] = 'ExclusiveHighCases'

    # Merge to bring in coefs/odds_ratio from featuresDf (same as HighCases path)
    featureImportance = featureImportance.merge(
        featuresDf, on=['feature', 'features_used'], how='left'
    )

    n_important = int(featureImportance['feature_importance'].sum())
    print(f'    {cohort} × {feature_str}: {n_successful}/{n_iterations} iterations OK  '
          f'| {n_important} features |z| > {z_threshold}')

    # ---- Plot z-score bar chart (reuse existing helper) ---------------------
    plot_and_save_top_features(
        z_scores, cohort,
        figPath, f'{feature_str}.ExclusiveHighCases'
    )

    # ---- Append to shared featureImportance CSV -----------------------------
    # Reindex to match the existing file's column order to prevent schema drift
    # when HighCases and ExclusiveHighCases have columns in different positions.
    if not os.path.exists(featureImportanceFile):
        with open(featureImportanceFile, mode='w', newline='') as fh:
            featureImportance.to_csv(fh, index=False)
    else:
        existing_cols = pd.read_csv(featureImportanceFile, nrows=0).columns.tolist()
        featureImportance = featureImportance.reindex(columns=existing_cols)
        with open(featureImportanceFile, mode='a', newline='') as fh:
            featureImportance.to_csv(fh, index=False, header=False)


def assign_groups(dfDict, featuresDf, trainingDataHigh, trainingDataLow,
                  high_thresh, low_thresh):
    groupDfList = []
    for g, df in dfDict.items():
        if 'IID' not in df.columns and 'Unnamed: 0' in df.columns:
            df.rename(columns={'Unnamed: 0': 'IID'}, inplace=True)
        df.set_index('IID', inplace=True)
        trainingDataHigh[g] = 0
        trainingDataHigh.loc[trainingDataHigh[f'bin_{g}'] > high_thresh, g] = 1
        trainingDataLow[g] = 0
        trainingDataLow.loc[trainingDataLow[f'bin_{g}'] < low_thresh, g] = 1
        df = df[featuresDf[featuresDf['model'] == g]['feature'].tolist()]
        groupDfList.append(df)
    return groupDfList


def main(phenoData, feature_scores_file, threshold, filter_strategy='tiered',
         majority_threshold=0.6, delta_threshold=0.5, use_all=False,
         run_exclusive=False, n_exclusive_iterations=1000,
         exclusive_z_threshold=1.0, filter_only=False, use_set='holdout'):

    scoresPath = f'{phenoData}/scores'
    figPath = f'{phenoData}/figures/importantCohortFeatures'

    featureImportanceFile = f"{scoresPath}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.csv"
    aucTableFile = f"{scoresPath}/importantCohortScores/performanceMetricsAcrossCohortsAndTrainingData.csv"
    filteredImportantFeaturesFile = f"{scoresPath}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv"

    if filter_only:
        print(f"\n[filter_only] Skipping training — reading existing results from:\n"
              f"  {featureImportanceFile}")
        if not os.path.exists(featureImportanceFile):
            raise FileNotFoundError(
                f"--filter_only requires an existing featureImportanceFile:\n"
                f"  {featureImportanceFile}"
            )
        # Jump directly to filtering section below
    else:
        # Reset output files so each full run starts clean rather than
        # appending to stale results from a previous run with a different schema.
        for _f in (featureImportanceFile, aucTableFile):
            if os.path.exists(_f):
                os.remove(_f)
        if use_set == 'holdout':
            combined_csv = f'{scoresPath}/combinedPRSGroups.holdout.filtered.csv'
        else:
            combined_csv = f'{scoresPath}/combinedPRSGroups.filtered.csv'
        print(f'[use_set={use_set}] Reading combined PRS from: {combined_csv}')
        combinedDf = pd.read_csv(combined_csv)
        models_bins_to_use = [m for m in combinedDf.columns if 'bin' in m]
        print('models to use ...', models_bins_to_use)

        models_to_use = [m.replace('bin_', '') for m in models_bins_to_use]
        print('models to use are : ', models_to_use)

        pattern = os.path.join(f'{scoresPath}/prsScores', "*mixed*")
        print('pattern to look for files is ..', pattern)

        tempFiles = glob.glob(pattern)
        if use_set == 'holdout':
            prsFiles = [f for f in tempFiles if 'holdout' in f]
        else:
            prsFiles = [f for f in tempFiles if 'holdout' not in f]
        prsFilesFiltered = [f for f in prsFiles if 'FromAll' not in f and not f.endswith('filtered.prs.csv')]
  
        models_sorted = sorted(models_to_use, key=len, reverse=True)
        filtered_files = []
        for filepath in prsFilesFiltered:
          filename = os.path.basename(filepath)
          print('[DEBUG] filename before second regex search : ',filename)
          for model in models_sorted:
            print('[DEBUG] model match in second regex search to account for main in epi+main : ',model)
            pat = rf'(^|[._\-\s]){re.escape(model)}([._\-\s]|$)'
            print('[DEBUG] regex pattern search : ',pat)
    #       if filename == f"{model}.mixed.prs.csv":
    #         filtered_files.append(filepath)
            if re.search(pat, filename):
              filtered_files.append(filepath)
              break
  
        dfDict = {}
        for file in filtered_files:
          print('[DEBUG] reading PRS from file ', file)
          data_type = file.split('.')[0].split('/')[-1]
          df = pd.read_csv(file)
          columns_to_drop = [col for col in df.columns if 'prs' in col]
          df.drop(columns=['PHENOTYPE'] + columns_to_drop, inplace=True)
          dfDict[data_type] = df
  
        featuresDf = pd.read_csv(feature_scores_file)
        featuresDf = featuresDf[
          (featuresDf['coefs'] != 0) &
          (~featuresDf['feature'].isin(['(Intercept)', 'SEX', 'age',
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5',
            'PC6', 'PC7', 'PC8', 'PC9', 'PC10']))
        ]
        featuresDf = featuresDf[featuresDf['model'].isin(models_to_use)]
  
        groupedDf = combinedDf[models_bins_to_use + ['IID', 'PHENOTYPE']].set_index('IID')

        high_thresh, low_thresh = detect_bin_thresholds(groupedDf, models_bins_to_use)

        high_risk_mask = pd.Series(False, index=groupedDf.index)
        low_risk_mask = pd.Series(False, index=groupedDf.index)
        for col in models_bins_to_use:
            high_risk_mask |= (groupedDf[col] > high_thresh)
            low_risk_mask |= (groupedDf[col] < low_thresh)

        trainingDataHigh = groupedDf[high_risk_mask & (groupedDf['PHENOTYPE'] == 2)]
        trainingDataLow = groupedDf[low_risk_mask & (groupedDf['PHENOTYPE'] == 1)]

        groups = dfDict.keys()
        groupDfList = assign_groups(dfDict, featuresDf, trainingDataHigh, trainingDataLow,
                                    high_thresh, low_thresh)
        featuresDf.rename(columns={'model': 'features_used'}, inplace=True)
  
        for trainingDataTuple in [(trainingDataHigh, 'HighCases'), (trainingDataLow, 'LowControls')]:
          trainingData = trainingDataTuple[0]
          trainingData_str = trainingDataTuple[1]
          probDfFile = f"{scoresPath}/importantCohortScores/probabilitiesAcrossCohortsAndFeaturesUsedInTraining.{trainingData_str}.csv"
  
          print(f'########################  TRAINING {trainingData_str}  ########################')
          probDf = pd.DataFrame()
  
          for data_type, df in dfDict.items():
            feature_str = data_type
            print(f'###################  FEATURE IMPORTANCE WITH {feature_str} FEATURES ###################')
  
            y = trainingData[groups]
            X = trainingData[['PHENOTYPE']].merge(df, left_index=True, right_index=True).drop(columns=['PHENOTYPE'])
            features = X.columns
  
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)
  
            model = MultiOutputClassifier(RandomForestClassifier(n_estimators=1000))
            model.fit(X_train, y_train)
  
            y_pred = model.predict(X_test)
            report_dict = classification_report(y_test, y_pred, target_names=groups, output_dict=True, zero_division=0)
            report_df = pd.DataFrame(report_dict).transpose().reset_index()
  
            for i, cl in enumerate(groups):
              proba_i = model.predict_proba(X_test)[i]
              proba_pos = proba_i[:, 1] if proba_i.shape[1] > 1 else proba_i[:, 0]
              probDfClass = pd.DataFrame(data=proba_pos)
              probDfClass.columns = [f'yHat.{cl}.features.{feature_str}']

              fpr, tpr, _ = roc_curve(y_test[cl], proba_pos)
              roc_auc = auc(fpr, tpr)
              report_df.loc[report_df['index'] == cl, 'auc'] = roc_auc
              print(f'AUC for {cl}: {roc_auc}')
  
              plt.figure()
              plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC for {cl} using {feature_str} features (area = {roc_auc:.2f})')
              plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
              plt.xlabel('FPR')
              plt.ylabel('TPR')
              plt.title(f'AUC for class {cl} using features {feature_str}')
              plt.savefig(f'{figPath}/{cl}.features.{feature_str}.{trainingData_str}.AUC.png')
              plt.close()
  
              explainer = fasttreeshap.TreeExplainer(model.estimators_[i], algorithm='auto', n_jobs=-1)
              shap_explainer = explainer(X_test, check_additivity=False)
              shap_values = shap_explainer.values
              # Binary RF: shape (n_samples, n_features, 2); single-class edge case: shape (..., 1)
              shap_class = shap_values[:, :, 1] if shap_values.shape[2] > 1 else shap_values[:, :, 0]

              shap_df = pd.DataFrame(data=shap_class, columns=features)
              create_summary_plots(figPath, shap_class, X_test, i, f'{cl}.{feature_str}.{trainingData_str}.')
  
              topFeatures, mean_shap_values, shap_z_scores = get_top_shap_features(shap_df, data_type, threshold)
              plot_and_save_top_features(shap_z_scores, i, figPath, data_type)

              featureImportance = pd.concat([shap_z_scores, mean_shap_values], axis=1)
              featureImportance.columns = ['z_score', 'mean_shap']
              featureImportance.dropna(inplace=True)
              featureImportance.reset_index(names="feature", inplace=True)
  
              featureImportance['feature_importance'] = 0
              featureImportance.loc[featureImportance['feature'].isin(topFeatures.index.tolist()), 'feature_importance'] = 1
  
              featureImportance['cohort'] = cl
              featureImportance['features_used'] = feature_str
              featureImportance['training_data_used'] = trainingData_str
              featureImportance = featureImportance.merge(featuresDf, on=['feature', 'features_used'], how='left')
  
              probDf = pd.concat([probDf, probDfClass], axis=1)
  
              if not os.path.exists(featureImportanceFile):
                with open(featureImportanceFile, mode='w', newline='') as f:
                  featureImportance.to_csv(f, index=False)
              else:
                with open(featureImportanceFile, mode='a', newline='') as f:
                  featureImportance.to_csv(f, index=False, header=False)
  
              report_df['features_used'] = feature_str
              report_df['training_data_used'] = trainingData_str
  
              if not os.path.exists(aucTableFile):
                with open(aucTableFile, mode='w', newline='') as f:
                  report_df.to_csv(f, index=False)
              else:
                with open(aucTableFile, mode='a', newline='') as f:
                  report_df.to_csv(f, index=False, header=False)
  
          probDf.index = X_test.index
          probDf.to_csv(probDfFile)

        # =========================================================================
        #  EXCLUSIVE HIGH-CASES: BOOTSTRAP GBM-SHAP PER COHORT
        #
        #  For each cohort find individuals exclusively high-risk there (bin > 8 in
        #  this cohort, NOT > 8 in any other), train n_exclusive_iterations
        #  HistGradientBoostedClassifier models each time drawing a fresh balanced
        #  control sample, extract SHAP values on cases, average across iterations,
        #  z-score, and append results to the shared featureImportanceFile under
        #  training_data_used = 'ExclusiveHighCases'.  These rows then participate
        #  in the same cross-cohort filtering step below as HighCases / LowControls.
        # =========================================================================
        if run_exclusive:
            print('\n########################  EXCLUSIVE HIGH-CASES BOOTSTRAP'
                  '  ########################')
            # Only use one bin col per loaded model — excludes combined, case_bin,
            # n_case_bins, prop_case_bins, and any column matching 'bin' via substring
            # (e.g. prs_comBINed). bin_combined is also excluded because it aggregates
            # all models and would make is_high_other True for nearly everyone.
            core_bin_cols = [f'bin_{m}' for m in dfDict.keys()
                             if f'bin_{m}' in groupedDf.columns]
            exclusive_groups = identify_exclusive_high_cases(groupedDf, core_bin_cols)

            for cohort, group_data in exclusive_groups.items():
                print(f'\n  --- Cohort: {cohort} ---')
                for data_type, feat_df in dfDict.items():
                    run_exclusive_high_bootstrap(
                        cohort               = cohort,
                        case_index           = group_data['case_index'],
                        control_index        = group_data['control_index'],
                        feature_df           = feat_df,
                        featuresDf           = featuresDf,
                        feature_str          = data_type,
                        n_iterations         = n_exclusive_iterations,
                        z_threshold          = exclusive_z_threshold,
                        figPath              = figPath,
                        featureImportanceFile= featureImportanceFile,
                        rng_seed             = 42,
                    )

    # =========================================================================
    #  CROSS-COHORT FILTERING with chosen strategy
    # =========================================================================
    # (reached either after training completes, or directly when --filter_only)
    importantTrainingDataFeatures = pd.read_csv(featureImportanceFile)
    importantTrainingDataFeatures.drop_duplicates(inplace=True)
    filteredResults = pd.DataFrame()

    for trainingData_str in importantTrainingDataFeatures['training_data_used'].unique():
        trainingDataFeatures = importantTrainingDataFeatures[
            importantTrainingDataFeatures['training_data_used'] == trainingData_str]
        feature_sets = trainingDataFeatures['features_used'].unique()

        # ExclusiveHighCases: each cohort was trained independently on its own
        # exclusive individuals. Two constraints for the filtered output:
        #   1. features_used == cohort: the feature set must be the one that
        #      actually defines this cohort's high-risk individuals; a feature
        #      important under another cohort's feature set is not useful downstream.
        #   2. max-z per feature: one authoritative row per feature, assigned to
        #      the cohort where its bootstrap signal was strongest.
        # Cross-cohort prevalence (how many cohorts flagged each feature) is stored
        # in n_exclusive_cohorts / exclusive_cohorts for biological interpretation.
        if trainingData_str == 'ExclusiveHighCases':
            important = trainingDataFeatures[
                trainingDataFeatures['feature_importance'] == 1
            ].copy()
            if important.empty:
                print(f"\n  ExclusiveHighCases: no features with feature_importance=1.")
                continue

            # Annotate cross-cohort prevalence from ALL important rows before filtering
            cohort_list = (important.groupby('feature')['cohort']
                           .agg(lambda s: ','.join(sorted(s.unique()))))
            n_cohorts_map = important.groupby('feature')['cohort'].nunique()

            # Keep only rows where the feature set matches the cohort
            matched = important[important['features_used'] == important['cohort']].copy()
            if matched.empty:
                print(f"\n  ExclusiveHighCases: no features with features_used == cohort.")
                continue

            # Among matching rows, keep the single row with max z_score per feature
            best_idx = matched.groupby('feature')['z_score'].idxmax()
            results_df = matched.loc[best_idx].copy()

            # Attach cross-cohort annotations
            results_df['n_exclusive_cohorts'] = results_df['feature'].map(n_cohorts_map)
            results_df['exclusive_cohorts'] = results_df['feature'].map(cohort_list)
            results_df['training_data_used'] = trainingData_str

            multi = (results_df['n_exclusive_cohorts'] > 1).sum()
            print(f"\n  ExclusiveHighCases: {len(important)} feature_importance=1 rows "
                  f"→ {len(matched)} after features_used==cohort "
                  f"→ {len(results_df)} after max-z assignment "
                  f"({multi} features important in >1 exclusive cohort).")
            filteredResults = pd.concat([filteredResults, results_df], ignore_index=True)
            continue

        # Print diagnostics so you can see where the sparsity originates
        print_feature_diagnostics(trainingDataFeatures, feature_sets, threshold=0.5)

        print(f"\nApplying filter strategy: '{filter_strategy}'  (training={trainingData_str})")

        if filter_strategy == 'strict':
            results_df = find_opposite_sign_features(trainingDataFeatures, feature_sets, threshold=0.5)

        elif filter_strategy == 'majority':
            results_df = find_majority_sign_features(
                trainingDataFeatures, feature_sets,
                threshold=0.5, majority_threshold=majority_threshold)

        elif filter_strategy == 'differential':
            results_df = find_differential_features(
                trainingDataFeatures, feature_sets,
                threshold=0.5, delta_threshold=delta_threshold)

        elif filter_strategy == 'tiered':
            results_df = find_tiered_features(
                trainingDataFeatures, feature_sets,
                threshold=0.5,
                majority_threshold=majority_threshold,
                delta_threshold=delta_threshold)
        else:
            raise ValueError(f"Unknown filter_strategy '{filter_strategy}'. "
                             "Choose from: strict, majority, differential, tiered")

        print(f"  → {len(results_df)} features survived ({trainingData_str})")

        results_df['training_data_used'] = trainingData_str
        if 'odds_ratio' in results_df.columns:
            results_df = results_df.sort_values('odds_ratio', ascending=False)

        filteredResults = pd.concat([filteredResults, results_df], ignore_index=True)

    filteredResults.to_csv(filteredImportantFeaturesFile, index=False)
    print(f"\nSaved filtered results → {filteredImportantFeaturesFile}")
    print(f"Total features across all training data: {len(filteredResults)}")

    # Summary table by strategy tier (if tiered mode)
    if filter_strategy == 'tiered' and 'specificity_tier' in filteredResults.columns:
        print("\nTier breakdown:")
        print(filteredResults.groupby(['training_data_used', 'specificity_tier']).size()
              .rename('n_features').reset_index().to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculating important features in statistically significant cohorts ....")
    parser.add_argument("--pheno_data", help="path to pheno results directory")
    parser.add_argument("--feature_scores_file_filtered", help="path to feature scores file")
    parser.add_argument("--threshold", type=float, default=1,
                        help="Z-score threshold for identifying important features per cohort (default: 1.99)")
    parser.add_argument("--filter_strategy", default='tiered',
                        choices=['strict', 'majority', 'differential', 'tiered'],
                        help=(
                            "Cross-cohort filtering strategy (default: tiered).\n"
                            "  strict       – original: positive in matching, ALL others negative\n"
                            "  majority     – positive in matching, negative in >= majority_threshold fraction of others\n"
                            "  differential – matching z is highest AND delta vs mean-others >= delta_threshold\n"
                            "  tiered       – runs all three and labels each feature by strictest tier it qualifies for"
                        ))
    parser.add_argument("--majority_threshold", type=float, default=0.6,
                        help="For majority strategy: fraction of other cohorts that must be negative (default: 0.6)")
    parser.add_argument("--delta_threshold", type=float, default=0.5,
                        help="For differential strategy: min z-score delta above mean of other cohorts (default: 0.5)")
    parser.add_argument("--run_exclusive", action='store_true', default=False,
                        help=(
                            "Run the exclusive high-cases bootstrap GBM-SHAP analysis. "
                            "Trains n_exclusive_iterations HistGradientBoosted models per "
                            "cohort × feature-set, each with a fresh balanced control "
                            "sample, averages SHAP z-scores across iterations, and appends "
                            "results to featureImportanceFile under "
                            "training_data_used='ExclusiveHighCases'. "
                            "Results participate in the same cross-cohort filtering step."
                        ))
    parser.add_argument("--n_exclusive_iterations", type=int, default=1000,
                        help=(
                            "Number of bootstrap iterations for the exclusive high-cases "
                            "GBM-SHAP analysis (default: 1000). Each iteration draws a "
                            "fresh balanced control sample, trains one GBM, and extracts "
                            "SHAP values on cases. 500 is a reasonable minimum; 1000 gives "
                            "stable z-score estimates in most cohort sizes."
                        ))
    parser.add_argument("--exclusive_z_threshold", type=float, default=1.0,
                        help=(
                            "SHAP z-score threshold for flagging feature_importance=1 in "
                            "the exclusive high-cases analysis (default: 1.0). Applied to "
                            "the z-score of the averaged per-feature mean SHAP across all "
                            "bootstrap iterations."
                        ))
    parser.add_argument("--filter_only", action='store_true', default=False,
                        help=(
                            "Skip all training and re-run only the cross-cohort filtering "
                            "step on the existing importantFeaturesAcrossCohortsAndTrainingData.csv. "
                            "Requires that file to already exist. Use this after manually "
                            "correcting that CSV without wanting to redo the full analysis."
                        ))

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")

    feature_scores_file = args.feature_scores_file_filtered or os.environ.get("FEATURE_SCORES_FILE_FILTERED")
    print(f"[PYTHON] Reading features from: {feature_scores_file}")

    threshold = os.environ.get("THRESHOLD")
    threshold = float(threshold) if threshold else args.threshold

#   pheno_data = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/combinedAnalysis'
#   feature_scores_file = f"{pheno_data}/scores/featureScoresReducedFinalModel.filtered.csv"

    if not pheno_data:
        raise ValueError("You must provide --pheno_data or set PHENO_DATA env var.")
    if not args.filter_only and not feature_scores_file:
        raise ValueError("You must provide --feature_scores_file_filtered or set FEATURE_SCORES_FILE_FILTERED env var.")

    print(f"Filter strategy        : {args.filter_strategy}")
    print(f"SHAP z threshold       : {threshold}")
    print(f"Majority threshold     : {args.majority_threshold}")
    print(f"Delta threshold        : {args.delta_threshold}")
    print(f"Run exclusive bootstrap: {args.run_exclusive}")
    print(f"Exclusive iterations   : {args.n_exclusive_iterations}")
    print(f"Exclusive z threshold  : {args.exclusive_z_threshold}")
    print(f"Filter only            : {args.filter_only}")

    main(
        pheno_data,
        feature_scores_file,
        threshold,
        filter_strategy=args.filter_strategy,
        majority_threshold=args.majority_threshold,
        delta_threshold=args.delta_threshold,
        run_exclusive=args.run_exclusive,
        n_exclusive_iterations=args.n_exclusive_iterations,
        exclusive_z_threshold=args.exclusive_z_threshold,
        filter_only=args.filter_only,
    )
  