#!/usr/bin/env python3

import pandas as pd
import numpy as np
#from helper.draw_plots import *

from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from scipy.stats import fisher_exact
import statsmodels.api as sm
import statsmodels.stats.contingency_tables as ct
from helper.draw_plots import *
from helper.download import *
from helper.data_wrangling import *
#from draw_optimized_auc_plots import *
from graph_grouped_prs_qq_plots import *

import glob


def standardize_prs(prsDf,prs_col):
	prsCopy = prsDf.copy()	
	scaler = StandardScaler().set_output(transform="pandas")
	scaled_model = scaler.fit(prsCopy[[prs_col]])
	scaled_feature = scaled_model.transform(prsCopy[[prs_col]])
	
	prsCopy[f'scaled_{prs_col}'] = scaled_feature
	return(prsCopy,scaled_model)


def calculate_prevalance(df):
	dfCopy = df.copy()
	dfCopy = dfCopy[dfCopy['PHENOTYPE'] == 2]
	dfCases = dfCopy.groupby(['bin']).count().rename(columns={'color':'total_cases'})[['total_cases']]
#	dfCases.rename(columns={'color':'total_cases'},inplace=True)[['total_cases']]
#	prevalenceDf = dfCases[['total_cases']]
	#get the count in each decile
	dfTotal = df.groupby(['bin']).count().rename(columns={'color':'total'})[['total']]
#	dfTotal.rename(columns={'color':'total'},inplace=True)[['total']]
	#get the prevalence in each decile
	dfTotal['prevalence'] = round((dfCases['total_cases'] / dfTotal['total']),2)
	dfTotal.reset_index(inplace=True)
	dfTotal = dfTotal[['prevalence','bin']]
	return(dfTotal)


def bin_prs(df,prs_column):
	color = df.iloc[0]['color']
	model = df.iloc[0]['model']
	
	dfCopy = df.copy()
	
	#sort scaled prs for entire dataset and break up into centiles
	dfCopy.sort_values([prs_column],inplace=True)
	
	#bin based on scaled iPRS (scaled_prs) into centiles
	dfCopy['bin'] = pd.qcut(dfCopy[prs_column], 10, labels=list(range(1,11)),duplicates='drop')
	dfCopy['bin_limits'] = pd.qcut(dfCopy[prs_column], 10,duplicates='drop')
	
	#round the bin limits to 1 decimal point
	dfCopy['bin_limits'] = dfCopy['bin_limits'].apply( lambda x : pd.Interval(left=round(x.left,1),right=round(x.right,1)))
	binLimits = dfCopy[['bin','bin_limits']].drop_duplicates()

	#prevalence df
	prevalenceDf = calculate_prevalance(dfCopy)	
	prevalenceDf = prevalenceDf.merge(binLimits,on=['bin'],how='left')
	prevalenceDf['bin'] = prevalenceDf['bin'].astype(float)
	prevalenceDf['color'] = color
	prevalenceDf['model'] = model
	dfCopy['bin'] = dfCopy['bin'].astype(float)
	return(prevalenceDf,dfCopy)


	
def apply_max_min(row,pheno_col,columns):
#	columns = [col for col in row.index if col != pheno_col]

	if row[pheno_col] == 2:  
		return row[columns].idxmax()  # Apply .idxmax() 
	else:
		return row[columns].idxmin()  # Apply .idxmin() to columns 
	
def assign_test_bin(row,prscr_col):
	cohort_bin = row[prscr_col]
	return(row[cohort_bin])


		
		
		
def main(pheno,scoresPath,figPath):
	
	data_str = 'Withcardio_main'
	#import the prs dataset and process based on which plots creating
	filePath = '/Users/kerimulterer/ukbiobank'
	figPath = f'{filePath}/{pheno}/tanigawaSet/figures/validation/sectionModels/paperFigures/filteredMainAllEpi'
	#columns to download in prs files
	columns=['prs','scaled_prs','PHENOTYPE','IID']
	
	#'PRScr_geno_with_cardio','PRScr_geno','PRScr_all'
	color_dict = {'main':'red','epi':'blue','epi+main':'purple','cardio':'#f5a142','PRScr_geno':'#f73bd7','PRScr_all':'#9c2287','all':'#c4771f'}
	marker_dict = {'main':'D','epi':'o','epi+main':'X','cardio':'v','PRScr_geno':'8','PRScr_all':'8','all':'v'}

	thresholdDf = pd.DataFrame()
		
		
		#instantiate files to download
	allFiles = glob.glob(f'{scoresPath}/*mixed.prs.csv')
	
	#cardioFile = [f for f in allFiles if 'Cardio' in f]
	
	#no covariate only data used in binning
	allFiles = [f for f in allFiles if 'covar' not in f]
	#remove the separate features from all model
	allFiles = [f for f in allFiles if 'FromAll' not in f]		
	
	
	outputPrevalenceFile = f'combinedPrevalencePRSGroups{data_str}.csv'
	outputCombinedGroupedFile = f'combinedPRSGroups{data_str}.csv'
	outputOddsRatioFile = f'combinedORPRSGroups{data_str}.csv'

	
	combinedPRSBinned = pd.DataFrame()
	prevalenceDf = pd.DataFrame()

	
	for file in allFiles:
		
		#main, epi, and epi+main are in the first position of each file name
#		data_type = file.split('/')[-1].split('.')[0]
		final_data_type = file.split('/')[-1].split('.')[0]
		
		#get the separateCardio with geno files
#		if 'Cardio' in file:
#			final_data_type = data_type+'+cardio'
#		else:
#			final_data_type = data_type
		
		#get the marker for the plots
		marker = marker_dict[final_data_type]

			
		try:	
			df = pd.read_csv(file,usecols=columns)
			df.set_index('IID',inplace=True)
		except ValueError: #the IID is an index so needs to be reindexed
			df = pd.read_csv(file,usecols=['prs','scaled_prs','PHENOTYPE','Unnamed: 0'])
			df.rename(columns={'Unnamed: 0':'IID'},inplace=True)
			df.set_index('IID',inplace=True)
		
			
		df['model'] = final_data_type
		df['color'] = color_dict[final_data_type]

		binned_prs = 'scaled_prs'
		prevalence, dfBinned = bin_prs(df,'scaled_prs')
		
		prevalence['marker'] = marker
		dfBinned.reset_index(inplace=True)# Rename the index to 'IID'
		prevalenceDf = pd.concat([prevalence,prevalenceDf],ignore_index=True)


		if combinedPRSBinned.empty: #create the df with the first datatype

			combinedPRSBinned = dfBinned[['IID','bin','prs',binned_prs,'PHENOTYPE']]
			combinedPRSBinned.rename(columns={'bin':f'bin_{final_data_type}'},inplace=True)
			combinedPRSBinned.rename(columns={f'{binned_prs}':f'{binned_prs}_{final_data_type}','prs':f'prs_{final_data_type}'},inplace=True)
			
		else:

			prsTemp = dfBinned[['IID','bin','prs',binned_prs,'PHENOTYPE']]
			prsTemp.rename(columns={'bin':f'bin_{final_data_type}'},inplace=True)
			prsTemp.rename(columns={binned_prs:f'{binned_prs}_{final_data_type}','prs':f'prs_{final_data_type}'},inplace=True)
			
			combinedPRSBinned = combinedPRSBinned.merge(prsTemp, on=['IID','PHENOTYPE'])
			
			
	####################### CALCULATE PRScr  ######################
	
	#### CALCULATE DIFFERENTLY WITH AND WITHOUT CARDIO

	bin_columns = ['bin_main','bin_epi','bin_epi+main','bin_cardio']	
	geno_columns = ['bin_main','bin_epi','bin_epi+main']
	
#	bin_columns = ['bin_main','bin_epi','bin_cardio']	
#	geno_columns = ['bin_main','bin_epi']
	
	
	combinedPRSBinned['PRScr_geno'] = combinedPRSBinned[['PHENOTYPE']+geno_columns].apply(apply_max_min, axis=1,args=('PHENOTYPE', geno_columns, ))
	combinedPRSBinned['PRScr_all'] = combinedPRSBinned[['PHENOTYPE']+bin_columns].apply(apply_max_min, axis=1,args=('PHENOTYPE', bin_columns, ))
	

	
		
	#calculated prevlance and bin for each PRScr cohort:

	for cr_cohort in ['PRScr_geno','PRScr_all']:
			
		marker = marker_dict[cr_cohort]
		
		combinedPRSBinned[f'bin'] = combinedPRSBinned.apply(assign_test_bin, axis=1, args=(cr_cohort,))

			
		#calculate prevalence for the bin 
		combinedPRSBinned['color'] = color_dict[cr_cohort]

		prevalence = calculate_prevalance(combinedPRSBinned)
		prevalence['marker'] = marker
		prevalence['model'] = cr_cohort
		prevalence['bin'] = prevalence['bin'].astype(float)
		prevalence['color'] = color_dict[cr_cohort]

		combinedPRSBinned.drop(columns=['color'],inplace=True)
		prevalenceDf = pd.concat([prevalenceDf,prevalence],ignore_index=True)
		
		combinedPRSBinned.rename(columns={'bin':f'bin_{cr_cohort}'},inplace=True)
			

	combinedPRSBinned.to_csv(f'{dataPath}/{outputCombinedGroupedFile}',index=False)
	prevalenceDf.to_csv(f'{dataPath}/{outputPrevalenceFile}',index=False)
	#calculate the raw PRS for each PRScr

	create_optimized_prevalence_plot(prevalenceDf,figPath,data_str)
	
	create_qq_plot_groups(combinedPRSBinned,f'{figPath}/correlationPlots')
			
	percentileOR = calculate_odds_ratio_for_prs(combinedPRSBinned,binned_prs)
	
	percentileOR.to_csv(f'{dataPath}/{outputOddsRatioFile}',index=False)
	

	
if __name__ == '__main__':
	
#   pheno = sys.argv[1]
	pheno = 'celiacDisease'
	
	main(pheno)

	