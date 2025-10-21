#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import statsmodels.stats.contingency_tables as ct
import statsmodels.api as sm
from scipy.stats import fisher_exact

def create_epi_df(epiDf,pairList,combo="sum"):
    '''input : epiDf with snps as columns + PHENOTYPE
            pairList : [snp1pair1,snp2pair1,snp1pair2,snp2pair2...snp1pairN,snp2pairN]'''
    epiArrayFinal = pd.DataFrame()
    
    for pair in pairList:
        if combo == "sum":  
            snps = epiDf[pair.split(',')].sum(axis=1)
        else:
            snps = epiDf[pair.split(',')].prod(axis=1)
        snps.columns = pair
        epiArrayFinal = pd.concat([epiArrayFinal,snps],axis=1)
    epiArrayFinal.columns = pairList 
    
    return(epiArrayFinal)


def get_epi_snps(epiFeatures):
    '''input : list[str:pair1,str:pair2..str:pairN]
    output : list : unique[str:snp1,str:snp2...str:snpN]'''
    epiSnps = []
    for pair in epiFeatures:
        epiSnps = epiSnps + pair.split(',')
    epiSnps = list(set(epiSnps))
    return(epiSnps)

def scale_data(df):
    # Initialize the StandardScaler
    scaler = StandardScaler()
    
    # Fit the scaler to your data (compute the mean and standard deviation)
    scaler.fit(df)
    
    # Transform the data using the fitted scaler
    scaled_data = scaler.transform(df)
    
    # Create a new DataFrame with the scaled data
    scaled_df = pd.DataFrame(scaled_data, columns=df.columns,index=df.index)
    
    return(scaled_df)

def compare_gene_env_to_genetic(scoresPath):
    inputFile = f'{scoresPath}/featureScoresReducedFinalModel.csv'
    df = pd.read_csv(inputFile)
    
    #columns = coefs,model,feature
    #iterate over all non-zero GxGxE features and get the first element in list to compare with G or GxG models
    gene_env = df[(df['model'] == 'cardio') & (df['coefs'] != 0) & (df['feature'].str.contains(','))]
    filterList = []
    for feature in gene_env['feature']:
        g = ','.join(feature.split(',')[1:])
        eDf = gene_env[gene_env['feature'] == feature]
        gDf = df[(df['feature'] == g) & (df['model'].isin(['epi','main']))]
        if gDf.empty:
            pass
        else:
            gValue = gDf.sort_values(['coefs'],ascending=False).head(1)['coefs'].values[0]
            eValue = eDf['coefs'].values[0]
            if gValue > 0: #for GxGxE features that are more predictive of risk
                if gValue > eValue:
                    filterList.append(feature)
            else:
                if eValue > gValue: #for GxGxE features that more predictive of controls
                    filterList.append(feature)
                
    dfFiltered = df[~df['feature'].isin(filterList)]
    inputFilePrefix = inputFile.split('.')[0]
    dfFiltered.to_csv(f'{inputFilePrefix}.filtered.csv')



def calculate_training_statistics(training_df):
    """
    Calculate mean and standard deviation for each PRS column in training data.
    
    Parameters:
    training_df: Training DataFrame with scaled PRS columns
    
    Returns:
    training_stats: Dictionary with mean and std for each PRS type
    """
    
    prs_columns = ['scaled_prs_main', 'scaled_prs_epi', 'scaled_prs_epi+main', 
                   'scaled_prs_cardio', 'scaled_prs_all']
    
    training_stats = {}
    
    print("Calculating training statistics for scaling...")
    
    for prs_col in prs_columns:
        if prs_col in training_df.columns:
            prs_values = training_df[prs_col].dropna()
            
            stats = {
                'mean': float(prs_values.mean()),
                'std': float(prs_values.std()),
                'count': len(prs_values),
                'min': float(prs_values.min()),
                'max': float(prs_values.max())
            }
            
            training_stats[prs_col] = stats
            
            print(f"  {prs_col}:")
            print(f"    Mean: {stats['mean']:.6f}")
            print(f"    Std: {stats['std']:.6f}")
            print(f"    Count: {stats['count']}")
        else:
            print(f"Warning: {prs_col} not found in training data")
            
    return training_stats

def scale_holdout_data_manually(holdout_df, training_stats):
    """
    Manually scale holdout data using training mean and standard deviation.
    
    Parameters:
    holdout_df: DataFrame with unscaled PRS columns
    training_stats: Dictionary with mean/std from calculate_training_statistics()
    
    Returns:
    holdout_scaled: DataFrame with scaled PRS columns added
    """
    
    # Base PRS column names (what we expect in holdout data)
    base_prs_columns = ['prs_main', 'prs_epi', 'prs_epi+main', 'prs_cardio', 'prs_all']
    
    holdout_scaled = holdout_df.copy()
    
    print("Manually scaling holdout data using training statistics...")
    
    for base_col in base_prs_columns:
        scaled_col = f'scaled_{base_col}'
        
        if base_col in holdout_df.columns and scaled_col in training_stats:
            
            # Get training statistics
            train_mean = training_stats[scaled_col]['mean']
            train_std = training_stats[scaled_col]['std']
            
            # Manual scaling: (x - mean) / std
            holdout_scaled[scaled_col] = (holdout_df[base_col] - train_mean) / train_std
            
            # Calculate statistics for verification
            scaled_values = holdout_scaled[scaled_col].dropna()
            
            print(f"  {base_col} -> {scaled_col}")
            print(f"    Training mean: {train_mean:.6f}, std: {train_std:.6f}")
            print(f"    Holdout scaled mean: {scaled_values.mean():.6f}, std: {scaled_values.std():.6f}")
            print(f"    Holdout scaled range: [{scaled_values.min():.3f}, {scaled_values.max():.3f}]")
            
        elif base_col not in holdout_df.columns:
            print(f"Warning: {base_col} not found in holdout data")
        elif scaled_col not in training_stats:
            print(f"Warning: No training statistics found for {scaled_col}")
            
    return holdout_scaled

def calculate_odds_centile_vs_all(df,decile_col,decile):
    '''calculate odds ratios using fisher exact test with pvalue and error bars
    
    input : df columns = individual group prs, PHENOTYPE,group decile
            decile for which OR is calculated
            prs_col is tuple (model_decile (float),PHENOTYPE (int),scaled_prs_model (float),model(str))
    
    output : odds ratio for decile with CI lower, CI upper
    
    a =  # High-risk group cases
    b =  # High-risk controls
    c =  # Rest of population cases
    d =  # Rest of population controls
    
    '''
    
    dfCopy = df.copy()
#	if decile == 25:
#		dfCopy['decile_group'] = (dfCopy[decile_col] <= decile).astype(int)
#	else:
    dfCopy['decile_group'] = (dfCopy[decile_col] > decile).astype(int)
    
    
    highRiskCases = dfCopy[(dfCopy['decile_group'] == 1) & (dfCopy['PHENOTYPE'] == 2)]
    highRiskControls = dfCopy[(dfCopy['decile_group'] == 1) & (dfCopy['PHENOTYPE'] == 1)]
    restCases = dfCopy[(dfCopy['decile_group'] == 0) & (dfCopy['PHENOTYPE'] == 2)]
    restControls = dfCopy[(dfCopy['decile_group'] == 0) & (dfCopy['PHENOTYPE'] == 1)]
    
    a = highRiskCases.shape[0]
    b = highRiskControls.shape[0]
    c = restCases.shape[0]
    d = restControls.shape[0]
    
    
    
    #get confidence intervals
    # Create a contingency table
    np_table = np.array([[a, b], [c, d]])
    
    sm_table = sm.stats.Table.from_data(np_table)
    
    # Calculate the odds ratio and confidence intervals
    table = ct.Table2x2(np_table)
    
    
    odds_ratio, p_value_fisher = fisher_exact(np_table)
    
    
    rslt = sm_table.test_nominal_association()
    p_value = rslt.pvalue
    if decile > 100:
        decile = int(decile / 10)
        
    if decile % 5 == 0:
        print('****************************************')
        print('')
        print('looking at odds ratio for the bottom ',decile, ' v top')
        print('number of high risk cases : ',a)
        print('number of high risk controls : ',b)
        print('number of low risk cases : ',c)
        print('number of low risk controls : ',d)
        
        print('pvalue for fishers exact = ',p_value_fisher)
        print('odds ratio for fishers exact = ',odds_ratio)
        print('pvalue for chi nominal association = ',p_value)
        
    # Perform chi-squared test
#	chi2_stat, p_value, dof, expected = ct.test_nominal_association(np_table)
#	print('pvalue for chi squared = ',p_value)
        
    odds_ratio_sm = table.oddsratio
    ci_low, ci_upp = table.oddsratio_confint()
    
    return(odds_ratio_sm,p_value_fisher,ci_low,ci_upp)

def calculate_odds_ratio_for_prs(df,binned_prs,prscr=False):
    '''calculate the odds ratio for each group based on the binned_prs 
    input : dataframe columns = ['IID', 'PHENOTYPE', 'prs_calculation', 'scaled_prs_calculation', 'color', 'bin_epi+main', 'bin_main', 'bin_epi',
        'prs_PRScr_model', 'scaled_prs_PRScr_model']
        
    output : dataframe columns = [percentile, odds_ratios, CI_upper, CI_lower,model]

        
        '''
    # run a logistic regression model using the PRS and phenotype
    # capture the beta coefficient
    # convert to odds ratios
    
    prs_columns = [col for col in df.columns if binned_prs in col]
#	prs_columns = [col for col in prs_columns if 'epi+main' not in col]
    cardio_prs = [col for col in prs_columns if 'cardio' in col]
    geno_prs = [col for col in prs_columns if 'cardio' not in col]
    
    
    percentileOR = pd.DataFrame()
    
    #data cleaning for columns
    for prs_col in prs_columns:
        print('')
        print('###############################################')
        print('data type = ',prs_col)
#		if prs_col == 'PRScr_mix':
#			bin_col = 'bin_PRScr_mix'
#		else:
        bin_col = prs_col.replace(binned_prs,'bin')
        try:
            df.sort_values(prs_col,inplace=True)
            # Assign numeric labels
            bins = pd.qcut(df[prs_col], 100, labels=list(range(1,101)),duplicates='drop')
        except IndexError:
            print('index error out of bounds when creating centiles')
            print('model which threw the error was :',prs_col)
            print('')
            return(percentileOR)
        percentileDf = df[[prs_col,'PHENOTYPE']]
        percentileDf[bin_col] = bins
#       percentileDf = percentileDf.rename(columns={prs_col:bin_col}).merge(df[['IID',prs_col,'PHENOTYPE']],left_index=True,right_index=True)
        
        ##################################################################################
        #                       CALCULATE TOP X VERSUS BOTTOM X                          #
        ##################################################################################
        
        
        percentileORTemp = pd.DataFrame()
        
        centiles = []
        centile_or = []
        centile_pvalues = []
        centile_ci_uppers = []
        centile_ci_lowers = []
        #compare the top and bottom centiles
        for centile in range(1,21):
            centiles.append(centile)
#			topBottomDf = percentileDf[(percentileDf[bin_col] <= centile) | (percentileDf[bin_col] > centile)]
            odds_ratio,p_value,ci_low,ci_upp = calculate_odds_centile_vs_all(percentileDf,bin_col,100-centile)
            centile_or.append(odds_ratio)
            centile_pvalues.append(p_value)
            centile_ci_uppers.append(ci_upp)
            centile_ci_lowers.append(ci_low)
            if centile % 5 == 0:
                
                print('centile # = ',centile)
                print('odds ratio = ',odds_ratio)
                
        percentileORTemp['percentiles'] = centiles
        percentileORTemp['odds_ratios'] = centile_or
        percentileORTemp['p_values'] = centile_pvalues
        percentileORTemp['CI_Upper'] = centile_ci_uppers
        percentileORTemp['CI_lower'] = centile_ci_lowers
        percentileORTemp['model'] = prs_col.replace(f'{binned_prs}_','')
        percentileOR = pd.concat([percentileORTemp,percentileOR],ignore_index=True)
        
    if prscr:
        for bin_col in ['bin_PRScr_all']:
#		for bin_col in ['bin_PRScr_all','bin_PRScr_geno']:
            print('')
            print('###############################################')
            print('data type = ',bin_col)
            
            percentileORTemp = pd.DataFrame()
            
            centiles = []
            centile_or = []
            centile_pvalues = []
            centile_ci_uppers = []
            centile_ci_lowers = []
#				for centile in [9,10]:
            for centile in range(1,21):
                centiles.append(centile)
#					topBottomDf = df[(df[bin_col] < centile*10) | (df[bin_col] >= centile*10)]
                odds_ratio,p_value,ci_low,ci_upp = calculate_odds_centile_vs_all(df,bin_col,1000-centile*10)
                
                
                if centile % 50 == 0:
                    print('centile # = ',centile)
                    print('odds ratio = ',odds_ratio)
                    
                centile_or.append(odds_ratio)
                centile_pvalues.append(p_value)
                centile_ci_uppers.append(ci_upp)
                centile_ci_lowers.append(ci_low)
                
            percentileORTemp['percentiles'] = centiles
            percentileORTemp['odds_ratios'] = centile_or
            percentileORTemp['p_values'] = centile_pvalues
            percentileORTemp['CI_Upper'] = centile_ci_uppers
            percentileORTemp['CI_lower'] = centile_ci_lowers
            percentileORTemp['model'] = bin_col.replace(f'bin_','')
            
            percentileOR = pd.concat([percentileORTemp,percentileOR],ignore_index=True)
        
        
    return(percentileOR)

#def rank_gene_env_features(geneEnvShapleyFile,threshold=2):
#   '''
#   input : df with columns [envGeneticFeature,shap_zscore,env_type,geneticFeature,envFeature,main_E,epistatic]
#   
#   oupt : df with gene_environment_features ranked with Shapley Z scores > 2
#       
#   '''
#   #download data
#   df = pd.read_csv(geneEnvShapleyFile)
#   
#   if 'Full' in geneEnvShapleyFile:
#       pass
#   else:
#       newFile = geneEnvShapleyFile.split('.')[0]
#       newFile = f'{newFile}Full.csv'
#       df.to_csv(newFile,index=False)
#   
#   #get the largest shap_zscore for duplicated values
#   df.sort_values(['envGeneticFeature','shap_zscore'],ascending=False,inplace=True)
#   
#   #drop duplicates in df
#   df.drop_duplicates(subset=['envGeneticFeature'],keep='first',inplace=True)
#   
#   #get the epistatic interactions
#   epiDf = df[df['epistatic'] == 1]
#   
#   #get the main effects
#   mainDf = df[df['main_E'] == 1]
#   mainDf.loc[mainDf['envFeature'].isna(),'envFeature'] = mainDf['envGeneticFeature']
#   
#   importantFeatures = epiDf[epiDf['shap_zscore'] > threshold]
#   
#   finalFeatures = pd.concat([importantFeatures,mainDf],ignore_index=True)
#   
#   finalFeatures.to_csv(geneEnvShapleyFile,index=False)
#   
#   return finalFeatures


#   
#   
#
#
#if __name__ == '__main__':
#   
#   main('/Users/kerimulterer/prsInteractive/results/type2Diabetes_test/epiFiles')


    