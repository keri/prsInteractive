#!/usr/bin/env python3
# Keri Multerer 2026
# Two-stage feature reduction: refactor of sklearnSectionModelsScoreTrain.py.
#
# V1 cost per section: ComplementNB + GridSearchCV over 4 learning rates x 3
# folds (12 XGBoost fits) + TreeExplainer SHAP over the full test set.
# V2 cost per section: ComplementNB + 1 early-stopped XGBoost fit + native
# gain importance, with SHAP restricted to a subsample of informative sections.
#
# Feature-space reduction (stage 1) is handled separately by
# univariate_screen.py, which filters the snplist/epi file this script reads.
#
# Output files and column layouts are identical to V1, so every downstream
# script (filter_features_LD.py, calculate_top_features_in_cohort.py, the glm
# modelling steps) works unchanged.

import pandas as pd
import numpy as np
import time
import xgboost as xgb
from sklearn.naive_bayes import ComplementNB
from sklearn.metrics import (balanced_accuracy_score, f1_score, roc_auc_score,
                             matthews_corrcoef, log_loss, jaccard_score, hamming_loss)
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn import set_config
from scipy import stats
import pickle
import csv
import warnings
import os
import argparse

set_config(transform_output="pandas")
warnings.simplefilter(action='ignore')

from helper.download import get_dataset, get_epi_columns, get_columns
from helper.data_wrangling import *


def convert_log_prob_to_odds(df):
    '''input df columns = ['log_prob_no_pheno','log_prob_yes_pheno']
       output : columns added ['odds_yes_pheno','odds_no_pheno']'''
    prob_no_pheno = np.exp(-(df['log_prob_no_pheno']))
    prob_yes_pheno = np.exp(-(df['log_prob_yes_pheno']))

    df['odds_yes_pheno'] = (df['log_prob_yes_pheno']/df['log_prob_no_pheno']) / (prob_yes_pheno/prob_no_pheno)
    df['odds_no_pheno'] = (df['log_prob_no_pheno']/df['log_prob_yes_pheno']) / (prob_no_pheno/prob_yes_pheno)

    df.sort_values(['odds_yes_pheno'], ascending=False, inplace=True)
    return(df)


def train_models(X, y, modelPath, pheno, data_type, i, learning_rate, n_estimators,
                 early_stopping_rounds, n_jobs, seed):
    '''input : X from training dataset and y
    output : fitted imputer, ComplementNB and XGBoost models, pickled to modelPath'''

    print('training models .....')
    print('shape of X in training dataset = ', X.shape)
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean').fit(X)
    Ximp = imp_mean.transform(X)
    pickle.dump(imp_mean, open(f'{modelPath}/imp_mean_{data_type}_{i}.pkl', 'wb'))

    ##########################
    # Naive Bayes
    ##########################

    st = time.time()
    clfNVB = ComplementNB().fit(Ximp, y)
    en = time.time()
    print('time if took to train naive bayes = ', (en-st)/60, ' minutes')

    ##########################
    # Gradient boosted classifier - single early-stopped fit
    ##########################
    # V1 ran GridSearchCV over 4 learning rates x 3 folds (12 fits per section).
    # Early stopping on a held-out slice of the training data tunes the number
    # of trees instead, which is the parameter that actually mattered, at 1/12
    # of the cost.
    st = time.time()
    Xtr, Xval, ytr, yval = train_test_split(
        X, y, test_size=0.2, random_state=seed,
        stratify=y if len(np.unique(y)) > 1 else None)

    clfHGB = xgb.XGBClassifier(
        tree_method='hist',
        max_bin=256,
        colsample_bytree=0.8,
        learning_rate=learning_rate,
        n_estimators=n_estimators,
        early_stopping_rounds=early_stopping_rounds,
        eval_metric='auc',
        n_jobs=n_jobs,
        random_state=seed,
    )
    clfHGB.fit(Xtr, ytr, eval_set=[(Xval, yval)], verbose=False)
    best_iter = getattr(clfHGB, 'best_iteration', None)
    en = time.time()
    print(f'time if took to train XGBoost = {(en-st)/60} minutes '
          f'(best_iteration={best_iter})')

    pickle.dump(clfNVB, open(f'{modelPath}/sklearnNaiveBayes_{data_type}_{i}.pkl', 'wb'))
    pickle.dump(clfHGB, open(f'{modelPath}/xgBoostClassifier_{data_type}_{i}.pkl', 'wb'))
    return(imp_mean, clfNVB, clfHGB)


def rank_features_by_gain(clfHGB, feature_names):
    '''Native XGBoost gain importance, standardised to a z-score so the
    threshold semantics match the SHAP z-scores V1 produced.'''
    booster = clfHGB.get_booster()
    gains = booster.get_score(importance_type='gain')
    # booster keys are the feature names when fitted from a DataFrame
    series = pd.Series({f: gains.get(f, 0.0) for f in feature_names})
    nonzero = series[series > 0]
    if len(nonzero) < 2:
        return pd.Series(dtype=float), series
    zscores = pd.Series(stats.zscore(nonzero), index=nonzero.index)
    return zscores, series


def calculate_shap_subsample(clfHGB, X, threshold, max_rows, seed):
    '''SHAP on a subsample. V1 explained every row of the test set; the
    feature ranking is a mean over rows, which converges long before the full
    set is used, so a capped subsample gives the same ordering far cheaper.
    Returns (topFeatures, zscores) or (None, None) if SHAP is unavailable.'''
    try:
        import shap
    except ImportError:
        print('shap not installed; falling back to gain importance only')
        return None, None

    Xs = X.sample(n=min(max_rows, len(X)), random_state=seed) if len(X) > max_rows else X
    try:
        explainer = shap.TreeExplainer(clfHGB, feature_perturbation='interventional',
                                       model_output='raw')
        shap_values = explainer.shap_values(Xs)
    except Exception as e:
        print(f'SHAP failed ({e}); falling back to gain importance only')
        return None, None

    if isinstance(shap_values, list):
        shap_values = shap_values[1]

    shap_df = pd.DataFrame(data=shap_values, columns=Xs.columns, index=Xs.index)
    shap_df = shap_df.loc[:, (shap_df != 0).any(axis=0)]
    if shap_df.shape[1] < 2:
        return None, None

    means = shap_df.mean().sort_values(ascending=False)
    zscores = pd.Series(stats.zscore(means), index=means.index)
    topFeatures = zscores[(zscores > threshold) | (zscores < -threshold)].sort_values()
    return topFeatures, zscores


def score_models(X, y, pheno, data_type, modelFile, i, imp_mean, clfNVB, clfHGB,
                 figPath, threshold, use_shap, shap_max_rows, seed):
    '''score both models on the test set and rank the section's features'''
    print('scoring models .....')
    st = time.time()

    Ximp = imp_mean.transform(X)

    #########################
    # Gradient Boosted
    #########################
    score = clfHGB.score(X, y)
    yHat = clfHGB.predict(X)
    proba = clfHGB.predict_proba(X)[:, 1]
    auc_hgb = roc_auc_score(y, proba)
    fields = ['xgboost', score, balanced_accuracy_score(y, yHat), auc_hgb,
              matthews_corrcoef(y, yHat), log_loss(y, proba),
              jaccard_score(y, yHat), hamming_loss(y, yHat),
              f1_score(y, yHat), data_type, i]
    with open(modelFile, mode='a') as f:
        csv.writer(f).writerow(fields)

    ##########################
    # Naive Bayes
    ##########################
    score = clfNVB.score(Ximp, y)
    yHat = clfNVB.predict(Ximp)
    proba_nb = clfNVB.predict_proba(Ximp)[:, 1]
    auc_nvb = roc_auc_score(y, proba_nb)
    fields = ['naive bayes', score, balanced_accuracy_score(y, yHat), auc_nvb,
              matthews_corrcoef(y, yHat), log_loss(y, proba_nb),
              jaccard_score(y, yHat), hamming_loss(y, yHat),
              f1_score(y, yHat), data_type, i]
    with open(modelFile, mode='a') as f:
        csv.writer(f).writerow(fields)

    rank_features = (auc_hgb > .51) or (auc_nvb > .51)
    print('time if took to score models = ', (time.time()-st)/60, ' minutes')

    dfSnps2 = pd.DataFrame()
    topFeatures = pd.DataFrame()

    if rank_features:
        dfSnps = pd.DataFrame()
        dfSnps['log_prob_no_pheno'] = clfNVB.feature_log_prob_[0, :]
        dfSnps['log_prob_yes_pheno'] = clfNVB.feature_log_prob_[1, :]
        dfSnps2 = convert_log_prob_to_odds(dfSnps)
        dfSnps2['feature'] = clfNVB.feature_names_in_

        gain_z, _ = rank_features_by_gain(clfHGB, list(X.columns))

        shap_top, shap_z = (None, None)
        if use_shap:
            shap_top, shap_z = calculate_shap_subsample(
                clfHGB, X, threshold, shap_max_rows, seed)

        if shap_z is not None:
            ranking, source = shap_z, 'shap'
            top = shap_top
        else:
            ranking, source = gain_z, 'gain'
            top = ranking[(ranking > threshold) | (ranking < -threshold)].sort_values()
        print(f'feature ranking source for section {i}: {source}')

        zdf = ranking.reset_index()
        zdf.columns = ['feature', 'shap_zscore']
        dfSnps2 = dfSnps2.merge(zdf, on='feature', how='left')

        if top is not None and len(top):
            topFeatures = top.rename('shap_zscore').to_frame()

    return(dfSnps2, topFeatures)


def main(pheno, withdrawalPath, pheno_path, training_path, test_path, epi_path,
         data_type, start, end, threshold, epi_combo, section_size,
         learning_rate, n_estimators, early_stopping_rounds, n_jobs,
         use_shap, shap_max_rows, seed):

    print('starting analysis .....')
    n = section_size

    scoresPath = f'{pheno_path}/scores'
    modelPath = f'{pheno_path}/models'
    figPath = f'{pheno_path}/figures'

    models = pd.DataFrame(columns=['model', 'test score', 'balanced score', 'auc',
                                   'matthews_corrcoef', 'log_loss', 'jaccard_score',
                                   'hamming_loss', 'f1_score', 'data_type', 'iteration'])
    modelFile = f'{scoresPath}/sklearnModelScoresSections.csv'
    if not os.path.exists(modelFile):
        print(f'creating model scores file : {modelFile} ...')
        with open(modelFile, mode='w', newline='') as f:
            models.to_csv(f, index=False)

    print('pheno files are being read from: ')
    print(pheno_path)

    full_columns = get_columns(training_path)
    start = int(start)
    end = int(end)

    featureScoresFile = f'{scoresPath}/featureScores.csv'
    importantFeaturesFile = f'{scoresPath}/importantFeaturesPostShap.csv'

    if not os.path.exists(featureScoresFile):
        allModelFeatures = pd.DataFrame(columns=['log_prob_no_pheno', 'log_prob_yes_pheno',
                                                 'odds_no_pheno', 'odds_yes_pheno',
                                                 'model#', 'feature', 'model', 'shap_zscore'])
        print(f'creating features scores file : {featureScoresFile} ...')
        with open(featureScoresFile, mode='w', newline='') as f:
            allModelFeatures.to_csv(f, index=False)

    if not os.path.exists(importantFeaturesFile):
        importantFeaturesShap = pd.DataFrame(columns=['feature', 'shap_zscore', 'data_type'])
        print(f'creating important features file : {importantFeaturesFile} ...')
        with open(importantFeaturesFile, mode='w', newline='') as f:
            importantFeaturesShap.to_csv(f, index=False)

    for i in range(start, end+1):
        istart = (i-1)*n
        if data_type == 'main':
            sectionSnps = full_columns[istart+6:istart+n]
            sectionPairs = [x for x in sectionSnps if pd.notna(x)]
            sectionSnps = sectionPairs
        else:
            print('epi_path to get filtered epi pairs ..', epi_path)
            epiColumns = get_epi_columns(epi_path)
            sectionPairs = epiColumns[istart:istart+n]
            sectionPairs = [x for x in sectionPairs if pd.notna(x)]
            sectionSnps = get_epi_snps(sectionPairs)

        print('number of snps in this section = ', len(sectionSnps))
        if len(sectionSnps) == 0:
            print(f'section {i} is empty, skipping')
            continue

        mainArray = get_dataset(training_path, withdrawalPath, sectionSnps, use_chunking=True)
        y = mainArray['PHENOTYPE'] - 1
        Xmain = mainArray.drop(columns=["PHENOTYPE"])

        if data_type != 'main':
            Xmain = create_epi_df(Xmain, sectionPairs, combo=epi_combo)

        print('Xmain array = ', Xmain.shape)
        print('training section model ....')

        if Xmain.shape[1] > 0:
            imp_mean, clfNVB, clfHGB = train_models(
                Xmain, y, modelPath, pheno, data_type, i, learning_rate,
                n_estimators, early_stopping_rounds, n_jobs, seed)

            testArray = get_dataset(test_path, withdrawalPath, sectionSnps, use_chunking=True)
            yTest = testArray["PHENOTYPE"] - 1
            Xtest = testArray.drop(columns=['PHENOTYPE'])

            if data_type != 'main':
                Xtest = create_epi_df(Xtest, sectionPairs, combo=epi_combo)

            print('Test array shape = ', Xtest.shape)

            allModelFeatures, topFeatures = score_models(
                Xtest, yTest, pheno, data_type, modelFile, i, imp_mean, clfNVB,
                clfHGB, figPath, threshold, use_shap, shap_max_rows, seed)

            if allModelFeatures.empty:
                continue

            allModelFeatures['model#'] = i
            allModelFeatures['model'] = data_type

            if not topFeatures.empty:
                topFeatures = topFeatures.reset_index()
                topFeatures.columns = ['feature', 'shap_zscore']
                topFeatures['data_type'] = data_type
                with open(importantFeaturesFile, mode='a', newline='') as f:
                    topFeatures.to_csv(f, index=False, header=False)

            with open(featureScoresFile, mode='a', newline='') as f:
                allModelFeatures.to_csv(f, index=False, header=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="running V2 two-stage models for wrapper...")
    parser.add_argument("--pheno_data", help="Path to the input pheno data directory")
    parser.add_argument("--training_file", help="data file of training data")
    parser.add_argument("--epi_file", help="epi file with epi pairs to analyze")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--data_type", help="data type to analyze")
    parser.add_argument("--start", help="start job")
    parser.add_argument("--end", help="end job")
    parser.add_argument("--pheno", help="Phenotype to analyze")
    parser.add_argument("--withdrawal_path", help="Genetic withdrawal path for IDs")
    parser.add_argument("--threshold", type=float, default=1,
                        help="Z-score threshold for identifying important features")
    parser.add_argument("--epi_combo", default='sum',
                        help='method to use for combining epi interactions (default: sum)')
    # V2 tuning
    parser.add_argument("--section_size", type=int, default=3000)
    parser.add_argument("--learning_rate", type=float, default=0.1)
    parser.add_argument("--n_estimators", type=int, default=500)
    parser.add_argument("--early_stopping_rounds", type=int, default=20)
    parser.add_argument("--n_jobs", type=int, default=8)
    parser.add_argument("--no_shap", action='store_true',
                        help='rank features by XGBoost gain only (fastest)')
    parser.add_argument("--shap_max_rows", type=int, default=2000,
                        help='cap on rows explained by SHAP per section')
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")

    pheno = args.pheno or os.environ.get("PHENO")
    print(f"[PYTHON] Phenotype : {pheno}")

    data_type = args.data_type or os.environ.get("DATA_TYPE")
    print(f"data type : {data_type}")

    training_path = args.training_file or os.environ.get("TRAINING_PATH")
    print(f"training file : {training_path}")

    test_path = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_path}")

    withdrawal_path = args.withdrawal_path or os.environ.get("WITHDRAWAL_PATH")
    print(f"reading withdrawals from file : {withdrawal_path}")

    epi_file = args.epi_file or os.environ.get("EPI_FILE")
    print(f"epi file being read : {epi_file}")

    start = args.start or os.environ.get("START")
    print(f"start : {start}")

    end = args.end or os.environ.get("END")
    print(f"end : {end}")

    threshold = os.environ.get("THRESHOLD")
    threshold = float(threshold) if threshold else args.threshold
    print(f"analyzing top features based on z score : {threshold}")

    epi_combo = args.epi_combo or os.environ.get("EPI_COMBO")
    print(f"method to use for combining epi interactions : {epi_combo}")

    use_shap = not args.no_shap
    print(f"SHAP ranking enabled : {use_shap} (max {args.shap_max_rows} rows per section)")

    if not epi_file:
        raise ValueError("You must provide a data type code via --epi_file or set the EPI_FILE environment variable.")
    if not pheno_data:
        raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
    if not data_type:
        raise ValueError("You must provide a data type code via --data_type or set the DATA_TYPE environment variable.")
    if not training_path:
        raise ValueError("You must provide a data type code via --training_path or set the TRAINING_PATH environment variable.")
    if not test_path:
        raise ValueError("You must provide a data type code via --test_path or set the TEST_PATH environment variable.")
    if not start:
        raise ValueError("You must provide a data type code via --start or set the START environment variable.")
    if not end:
        raise ValueError("You must provide a data type code via --end or set the END environment variable.")
    if not epi_combo:
        raise ValueError("You must provide a epi_combo via --epi_combo or set the EPI_COMBO environment variable.")

    for sub in ("models", "scores", "figures"):
        os.makedirs(f"{pheno_data}/{sub}", exist_ok=True)

    main(pheno, withdrawal_path, pheno_data, training_path, test_path, epi_file,
         data_type, start, end, threshold, epi_combo, args.section_size,
         args.learning_rate, args.n_estimators, args.early_stopping_rounds,
         args.n_jobs, use_shap, args.shap_max_rows, args.seed)
