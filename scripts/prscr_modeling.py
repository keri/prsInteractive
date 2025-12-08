#!/usr/bin/env python3

#Keri Multerer 2024
#Modeling several PRS to obtain PRS cr score

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score, matthews_corrcoef, log_loss, jaccard_score, hamming_loss
from sklearn.model_selection import GridSearchCV
from sklearn.utils import resample
from multiprocessing import Pool

from scipy.stats import norm
from helper.data_wrangling import *
from helper.calculate_prs import *
import warnings

def bootstrap_lasso(X, y, model):
    X_boot, y_boot = resample(X, y)
    lasso_boot = model.fit(X_boot,y_boot)
    #   lasso_boot.fit(X_boot, y_boot)
    return (lasso_boot.coef_[0])

def calculate_confidence_interval(bootstrap_coefs,alpha=.05):
    # Calculate 95% confidence intervals
    lower_bound = np.percentile(bootstrap_coefs, 100 * (alpha / 2), axis=0)
    upper_bound = np.percentile(bootstrap_coefs, 100 * (1 - alpha / 2), axis=0)
    return(lower_bound,upper_bound)

def calculate_r_squared(y_pred,y_true):
    # Calculate Cox & Snell R^2
    epsilon=1e-10
    ll_model = -log_loss(y_true, y_pred, normalize=False)
    
    # Calculate the log-likelihood of the null model
    p_null = np.mean(y_true)  # Probability for the null model
    y_null = np.full_like(y_true, p_null)
    ll_null = -log_loss(y_true, y_null, normalize=False)
    
    # Maximum possible log-likelihood
    ll_max = -log_loss(y_true, y_true, normalize=False)  # Log-likelihood when model perfectly fits
    
    # Add a small constant to ll_max to prevent division by zero
    if ll_max == 0:
        ll_max += epsilon
        
    # Cox & Snell R-squared
    r2_cs = 1 - (ll_model / ll_null)
    
    # Nagelkerke R-squared
    try:
        r2_nagelkerke = r2_cs / (1 - (ll_null / ll_max))
        return r2_nagelkerke
    
    except ZeroDivisionError:
        print('division error, returning cox snell r squared')
        return r2_cs
    


def lasso_p_values_parallel(X, y, fitted_model, n_bootstrap=1000, n_jobs=18):
    """
    Calculate p-values for coefficients obtained from a Lasso model using bootstrapping in parallel.

    Args:
    - X (array-like): Features matrix.
    - y (array-like): Target vector.
    - alpha (float): Significance level (default: 0.05).
    - n_bootstrap (int): Number of bootstrap samples (default: 1000).
    - n_jobs (int): Number of parallel jobs. Set to -1 to use all available cores (default).

    Returns:
    - p_values (array): P-values for each coefficient.
    - ci : list of upper and lower confidence intervals
    - se (array): for each beta
    """

    coefs = fitted_model.coef_[0]

    # Bootstrap sampling and coefficient estimation using multiprocessing
    with Pool(n_jobs) as pool:
        #coefs_boot is a list arrays len(n_jobs) containing arrays of shape(n features) with a coefficient for each
        coefs_boot = pool.starmap(bootstrap_lasso, [(X, y, fitted_model) for _ in range(n_bootstrap)])
        
    lower_bound,upper_bound = calculate_confidence_interval(coefs_boot)
    # Calculate p-values
    # Calculate standard errors
    se = np.std(coefs_boot, axis=0)
    
    # Calculate t-statistics
    t_stats = coefs / se
    
    # Calculate p-values
    p_values = 2 * (1 - norm.cdf(np.abs(t_stats)))
    
    return(coefs,p_values,se,lower_bound,upper_bound)

def train_model(X_train,y_train,X_test,y_test):
    
    ############  train models #########
    parameters_lasso = [{'max_iter':[500,1000,2000]}]
    clfLasso = LogisticRegression(fit_intercept=False,random_state=24,solver='lbfgs')
    grid_search_lasso = GridSearchCV(estimator=clfLasso,param_grid=parameters_lasso,scoring='roc_auc',cv=5,n_jobs=18)
    grid_search_lasso.fit(X_train,y_train)
    
    coefs,p_values,se,lower_bound,upper_bound = lasso_p_values_parallel(X_train, y_train, grid_search_lasso.best_estimator_)
    
    w = np.array(grid_search_lasso.best_estimator_.coef_[0]).transpose()
    
    ####### test model #####
    
    score = grid_search_lasso.score(X_test,y_test)
    yHat = grid_search_lasso.predict(X_test)
    yProba = grid_search_lasso.predict_proba(X_test)[:, 1]
    balanced_score = balanced_accuracy_score(y_test,yHat)
    auc = roc_auc_score(y_test, yProba)
    mcc = matthews_corrcoef(y_test,yHat)
    logloss = log_loss(y_test,yProba)
    jscore = jaccard_score(y_test,yHat)
    hloss = hamming_loss(y_test,yHat)
    r_squared = calculate_r_squared(yProba,y_test)
    fields=[score,balanced_score,auc,mcc,logloss,jscore,hloss,r_squared]
    
    coefs = pd.DataFrame({'prs_cohort':X_train.columns,'beta':coefs})
    
    return(coefs,p_values,se,lower_bound,upper_bound,fields)

def calculate_prscr_mix(df,holdoutDf,figPath,scoresPath,use_epi_main=False,use_all=False):
    
    '''run mixed prs model using scaled prs from validation set. Use beta coefficients from models to calculate prscr mix in holdout set
    input: validation dataframe columns : [IID,PHENOTYPE,scaled_prs_main,scaled_prs_epi,scaled_prs_cardio,PC1,..PC10,SEX,age]
           
            holdout dataframe : columns = [IID,PHENOTYPE] + [scaled_prs from validation columns, PC1, .. PC10, SEX, age]

    output : dataframe : [PRScr_mix,prs_PRScr_mix,scaled_prs_PRScr_mix, bin_PRScr_mix,decile_bin_PRScr_mix]
                         index = IID
    
    '''
    dfCopy = df.copy()
    holdoutDfCopy = holdoutDf.copy()
    
    prs_cols = [col for col in dfCopy.columns if 'scaled_prs' in col]
    
    if not use_all:
        #remove the separate features from all model
        prs_cols = [col for col in prs_cols if 'all' not in col]
        
    if not use_epi_main:
        allFilesTemp = [f for f in prs_cols if 'epi+main' not in f]
    
  
    validation_prs_to_model = dfCopy[prs_cols]
    validation_prs_to_model = dfCopy.drop(columns=['PHENOTYPE'])
    
    print('prs columns to model : {prs_cols}')
    
    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(validation_prs_to_model, dfCopy[['PHENOTYPE']], test_size=0.3, random_state=42)
    
    
    coefs_temp,p_values,se,lower_bound,upper_bound,fields = train_model(X_train,y_train.values.ravel(),X_test,y_test.values.ravel())
    #get the beta values
    coefs = coefs_temp.copy()
    
    print('p_values :',p_values)
    print('se :',se)
    print('fields :',fields)
    
    coefs_temp['p_values'] = p_values
    coefs_temp['upper_CI'] = upper_bound
    coefs_temp['lower_CI'] = lower_bound
    coefs_temp['se'] = se
    coefs_temp['lower_CI'] = lower_bound
    
    #fields = [score,balanced_score,auc,mcc,logloss,jscore,hloss,r_squared]
    coefs_temp['score'] = fields[0]
    coefs_temp['balanced_score'] = fields[1]
    coefs_temp['auc'] = fields[2]
    coefs_temp['nagelkerke_rsquared'] = fields[-1]
    
    coefs_temp.to_csv(f'{scoresPath}/prscr_mix_scores.csv',index=False)
    ###################  calculate PRScr  ###################

#   coefs_temp = pd.read_csv(f'{scoresPath}/prscr_mix_scores.csv')
    
    
    
    #get the multipliers for the scaled_prs
    # Ensure the column names match the index in betas for alignment
    betas = coefs_temp[coefs_temp['prs_cohort'].str.contains('scaled_prs')]
    betas.set_index('prs_cohort',inplace=True)
    beta_multiplier = betas['beta']
    
    # Multiply each column in df_values by the corresponding multiplier
    prscr_holdout_df = holdoutDfCopy[prs_cols].multiply(beta_multiplier, axis=1)
    

    holdoutDfCopy['prs_cr_mix'] = prscr_holdout_df.sum(axis=1)
    
    #scaled PRS calculations
    scaled_prscr_mix = scale_data(holdoutDfCopy[['prs_cr_mix']])
    scaled_prscr_mix.columns = ['scaled_prs_cr_mix']

    
    #create calucate bins and use for prevalence and prs plots 
    scaled_prscr_mix.sort_values(['scaled_prs_cr_mix'],inplace=True)
    scaled_prscr_mix['decile_bin_cr_mix'] = pd.qcut(scaled_prscr_mix['scaled_prs_cr_mix'], 10, labels=list(range(1,11)),duplicates='drop')
    scaled_prscr_mix['bin_cr_mix'] = pd.qcut(scaled_prscr_mix['scaled_prs_cr_mix'], 1000, labels=list(range(1,1001)),duplicates='drop')
    
    
    scaled_holdout = holdoutDfCopy[['PHENOTYPE','prs_cr_mix']].merge(scaled_prscr_mix,left_index=True,right_index=True)
    scaled_holdout.rename(columns={'scaled_prs_cr_mix':'scaled_prs'},inplace=True)
    
    #calculate plots for prscr_mix
    create_prs_plots(scaled_holdout[['PHENOTYPE','scaled_prs']],'prscr_mix',figPath,scoresPath,'mixed')
    
    scaled_holdout.rename(columns={'scaled_prs':'scaled_prs_cr_mix'},inplace=True)

    
    return(scaled_holdout.drop(columns=['PHENOTYPE']))
    
    
def main(dataPath,cov_file):
    '''input : phenotype
    
    downloads the scaled prs for validation set to calculate beta values to be used in prs cr

    validation df: columns = [IID,prs_main, scaled_prs_main, prs_epi, scaled_prs_epi, prs_cardio, scaled_prs_cardio, prs_epi+main, scaled_prs_epi+main, PHENOTYPE]

    output: 
    prs_beta_df : beta for prs cohorts, pvalues, CI upper, CI lower, se
    model scores list : [score,balanced_score,auc,mcc,logloss,jscore,hloss,r_squared]
    '''
    figPath = f'{dataPath}/figures'
    scoresPath = f'{dataPath}/scores'
        
    covDf = pd.read_csv(cov_file,usecols=['IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','age','SEX'])
    covDf.set_index('IID',inplace=True)
    
    df = pd.read_csv(f'{scoresPath}/combinedPRSGroups.csv')
    df.set_index('IID',inplace=True)
    df = df.merge(covDf,right_index=True,left_index=True,how='left')
    
    holdout_df = pd.read_csv(f'{scoresPath}/combinedPRSGroups.holdout.csv')
    holdout_df.set_index('IID',inplace=True)
    
    ############### CHECK EPI + MAIN STATISTICS ##################
    statsTable = f"{pheno_data}/scores/prs_pairwise_comparisons.csv"
    
    ############ IS EPI+MAIN INCLUDED IN IMPORTANT FEATURE ANALYSIS  ##############
    stats = pd.read_csv(statsTable)
    
    epiMainStats = stats[(stats['model2'] == 'scaled_prs_epi+main') & (stats['delong_p_value'] > .05)]
    
    if len(epiMainStats) > 0:
        print('NOT USING EPI+MAIN IN ANALYSIS')
        print(epiMainStats)
        use_epi_main = False
    else:
        use_epi_main = True
    
    prsHoldout = calculate_prscr_mix(df,holdout_df,figPath,scoresPath,use_epi_main=use_epi_main,use_all=False)
    holdout_df = holdout_df.merge(prsHoldout,left_index=True,right_index=True,how='left')
    holdout_df.reset_index(inplace=True)
    
    holdout_df.to_csv(f'{scoresPath}/combinedPRSGroups.holdout.csv',index=False)
    
#       prsHoldout.to_csv(f'{trainingPath}/prs/reducedSHAP/holdout/combinedPRSGroupsWithcardio_main.holdout.csv')
    ###############  CREATE PRS COMBINED DATA WITH THIS METHOD ###############
    #columns =  bin_main

    
if __name__ == '__main__':
    __spec__ = None
    parser = argparse.ArgumentParser(description="Calculating prscr_mix ....")
    parser.add_argument("--pheno_data",help="path to pheno results directory")
    parser.add_argument("--covar_file",help="path to cleaned covariate file")
    
    
    args = parser.parse_args()
    
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
    covar_file = args.covar_file or os.environ.get("COVAR_FILE")
    print(f"[PYTHON] Reading covariate data from: {covar_file}")
    
    pheno = 'celiacDisease'
    pheno_data = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi'
    covar_file = '/Users/kerimulterer/prsInteractive/results/covar.csv'
    
    if not pheno_data:
        raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
        
    if not covar_file:
        raise ValueError("You must provide a covar data file via --covar_file or set the COVAR_FILE environment variable.")
        


    main(pheno_data,covar_file)

    