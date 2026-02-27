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
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
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

def assign_groups(dfDict, featuresDf, trainingDataHigh, trainingDataLow):
    groupDfList = []
    for g, df in dfDict.items():
        try:
            df.rename(columns={'Unnamed: 0': 'IID'}, inplace=True)
        except KeyError:
            pass
        df.set_index('IID', inplace=True)
        trainingDataHigh[g] = 0
        trainingDataHigh.loc[trainingDataHigh[f'bin_{g}'] > 8, g] = 1
        trainingDataLow[g] = 0
        trainingDataLow.loc[trainingDataLow[f'bin_{g}'] < 3, g] = 1
        df = df[featuresDf[featuresDf['model'] == g]['feature'].tolist()]
        groupDfList.append(df)
    return groupDfList


def main(phenoData, feature_scores_file, threshold, filter_strategy='tiered',
         majority_threshold=0.6, delta_threshold=0.5, use_all=False):

    scoresPath = f'{phenoData}/scores'
    figPath = f'{phenoData}/figures/importantCohortFeatures'

    featureImportanceFile = f"{scoresPath}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.csv"
    aucTableFile = f"{scoresPath}/importantCohortScores/performanceMetricsAcrossCohortsAndTrainingData.csv"
    filteredImportantFeaturesFile = f"{scoresPath}/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.Filtered.{filter_strategy}.csv"

    combinedDf = pd.read_csv(f'{scoresPath}/combinedPRSGroups.filtered.csv')
    models_bins_to_use = [m for m in combinedDf.columns if 'bin' in m]
    print('models to use ...', models_bins_to_use)
  
    models_to_use = [m.replace('bin_', '') for m in models_bins_to_use]
    print('models to use are : ', models_to_use)
  
    pattern = os.path.join(f'{scoresPath}/prsScores', "*mixed*")
    print('pattern to look for files is ..', pattern)
  
    tempFiles = glob.glob(pattern)
    prsFiles = [f for f in tempFiles if 'holdout' not in f]
    prsFilesFiltered = [f for f in prsFiles if 'FromAll' not in f]
  
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
  
    high_risk_mask = pd.Series(False, index=groupedDf.index)
    low_risk_mask = pd.Series(False, index=groupedDf.index)
    for col in models_bins_to_use:
      high_risk_mask |= (groupedDf[col] > 8)
      low_risk_mask |= (groupedDf[col] < 3)
  
    trainingDataHigh = groupedDf[high_risk_mask & (groupedDf['PHENOTYPE'] == 2)]
    trainingDataLow = groupedDf[low_risk_mask & (groupedDf['PHENOTYPE'] == 1)]
  
    groups = dfDict.keys()
    groupDfList = assign_groups(dfDict, featuresDf, trainingDataHigh, trainingDataLow)
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
        report_dict = classification_report(y_test, y_pred, target_names=groups, output_dict=True)
        report_df = pd.DataFrame(report_dict).transpose().reset_index()
  
        for i, cl in enumerate(groups):
          probDfClass = pd.DataFrame(data=model.predict_proba(X_test)[i][:, 1])
          probDfClass.columns = [f'yHat.{cl}.features.{feature_str}']
  
          fpr, tpr, _ = roc_curve(y_test[cl], model.predict_proba(X_test)[i][:, 1])
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
  
          shap_df = pd.DataFrame(data=shap_values[:, :, 1], columns=features)
          create_summary_plots(figPath, shap_values[:, :, 1], X_test, i, f'{cl}.{feature_str}.{trainingData_str}.')
  
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
    #  CROSS-COHORT FILTERING with chosen strategy
    # =========================================================================
    importantTrainingDataFeatures = pd.read_csv(featureImportanceFile)
    importantTrainingDataFeatures.drop_duplicates(inplace=True)
    filteredResults = pd.DataFrame()

    for trainingData_str in importantTrainingDataFeatures['training_data_used'].unique():
        trainingDataFeatures = importantTrainingDataFeatures[
            importantTrainingDataFeatures['training_data_used'] == trainingData_str]
        feature_sets = trainingDataFeatures['features_used'].unique()

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
    if not feature_scores_file:
        raise ValueError("You must provide --feature_scores_file_filtered or set FEATURE_SCORES_FILE_FILTERED env var.")

    print(f"Filter strategy   : {args.filter_strategy}")
    print(f"SHAP z threshold  : {threshold}")
    print(f"Majority threshold: {args.majority_threshold}")
    print(f"Delta threshold   : {args.delta_threshold}")

    main(
        pheno_data,
        feature_scores_file,
        threshold,
        filter_strategy=args.filter_strategy,
        majority_threshold=args.majority_threshold,
        delta_threshold=args.delta_threshold
    )
  