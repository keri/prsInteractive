#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from helper.draw_plots import *
from helper.download import *
from helper.data_wrangling import *
import glob
import os
import argparse

def combine_prs_files(allFiles):
	'''
	combined prs from all datatypes into one dataframe

	input: list of all files in directory with prs of interest
	output: dataframe with columns [prs_{type1} ..etc, PHENOTYPE, IID]

	'''
	prsList = []
	combinedPRS = pd.DataFrame()
	for file in allFiles:
		#main, epi, and epi+main are in the first position of each file name
#		data_type = file.split('/')[-1].split('.')[0]
		final_data_type = file.split('/')[-1].split('.')[0]
		prsList.append(final_data_type)
		try:	
			df = pd.read_csv(file,usecols=['prs','PHENOTYPE','IID'])
			df.set_index('IID',inplace=True)
		except ValueError: #the IID is an index so needs to be reindexed
			df = pd.read_csv(file,usecols=['prs','PHENOTYPE','Unnamed: 0'])
			df.rename(columns={'Unnamed: 0':'IID'},inplace=True)
			df.set_index('IID',inplace=True)
			
		df.rename(columns={'prs':f'prs_{final_data_type}'},inplace=True)
		if combinedPRS.empty:
			combinedPRS = df.copy()
		else:
			combinedPRS = combinedPRS.merge(df, on=['IID','PHENOTYPE'])
	return combinedPRS, prsList
			


def calculate_prevalance(df):
	dfCopy = df.copy()
	dfCopy = dfCopy[dfCopy['PHENOTYPE'] == 2]
	dfCases = dfCopy.groupby(['bin']).count().rename(columns={'model':'total_cases'})[['total_cases']]
#	dfCases.rename(columns={'color':'total_cases'},inplace=True)[['total_cases']]
#	prevalenceDf = dfCases[['total_cases']]
	#get the count in each decile
	dfTotal = df.groupby(['bin']).count().rename(columns={'model':'total'})[['total']]
#	dfTotal.rename(columns={'color':'total'},inplace=True)[['total']]
	#get the prevalence in each decile
	dfTotal['prevalence'] = round((dfCases['total_cases'] / dfTotal['total']),2)
	dfTotal.reset_index(inplace=True)
	dfTotal = dfTotal[['prevalence','bin']]
	return(dfTotal)

def assign_risk_bins(df, prs_column, n_bins=10):
		"""
		Assign risk bins for PRS scores with proper error handling.
		"""
	
		print(f"Creating {n_bins} risk bins for {prs_column}")
	
		# Check if we have enough unique values
		unique_vals = df[prs_column].nunique()
		print(f"Unique values in {prs_column}: {unique_vals}")
	
		if unique_vals < n_bins:
				print(f"Warning: Only {unique_vals} unique values, reducing bins to {unique_vals}")
				n_bins = unique_vals
			
		if n_bins <= 1:
				print("Cannot create meaningful bins - assigning all to bin 1")
				return pd.Series([1] * len(df), index=df.index, name='bin')
	
		# Create bins
		try:
				# First try with qcut
				bins = pd.qcut(df[prs_column], n_bins, duplicates='drop')
				actual_bins = len(bins.cat.categories)
			
				if actual_bins != n_bins:
						print(f"Note: {n_bins} bins requested, {actual_bins} bins created due to duplicate values")
					
				# Assign numeric labels
				bins = pd.qcut(df[prs_column], n_bins, 
											labels=list(range(1, actual_bins + 1)), 
											duplicates='drop')
			
				return bins
	
		except Exception as e:
				print(f"qcut failed: {e}")
				# Fallback to rank-based binning
				ranks = df[prs_column].rank(method='first')
				bins = pd.cut(ranks, n_bins, labels=list(range(1, n_bins + 1)))
				return bins
		


def bin_prs(df,prs_column,prs):

	
	dfCopy = df.copy()
	
	#sort scaled prs for entire dataset and break up into centiles
	dfCopy.sort_values([prs_column],inplace=True)
	
	#bin based on scaled iPRS (scaled_prs) into centiles
	dfCopy['bin'] = assign_risk_bins(dfCopy, prs_column)
	dfCopy['bin_limits'] = pd.qcut(dfCopy[prs_column], 10,duplicates='drop')
	binLimits = dfCopy[['bin','bin_limits']].drop_duplicates()
	dfCopy['model'] = prs
	
	#prevalence df
	prevalenceDf = calculate_prevalance(dfCopy)
	#round the bin limits to 1 decimal point
	prevalenceDf = prevalenceDf.merge(binLimits,on=['bin'],how='left')
	prevalenceDf['bin'] = prevalenceDf['bin'].astype(float)
	
	
	dfCopy['bin'] = dfCopy['bin'].astype(float)
	dfCopy.rename(columns={'bin':f'bin_{prs}'},inplace=True)
	
	prevalenceDf[f'bin_limits'] = prevalenceDf['bin_limits'].apply( lambda x : pd.Interval(left=round(x.left,1),right=round(x.right,1)))

	return(prevalenceDf,dfCopy[[f'bin_{prs}']])

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


		
def main(phenoPath):
	
	#import the prs dataset and process based on which plots creating
	scoresPath = f'{phenoPath}/scores'
	figPath = f'{phenoPath}/figures'
	
	#'PRScr_geno_with_cardio','PRScr_geno','PRScr_all'
#	color_dict = {'main':'red','epi':'blue','epi+main':'purple','cardio':'#f5a142','PRScr_geno':'#f73bd7','PRScr_all':'#9c2287','all':'#c4771f'}
	color_dict = {'main':'red','epi':'blue','epi+main':'purple','cardio':'#f5a142','all':'#c4771f'}
	marker_dict = {'main':'D','epi':'o','epi+main':'X','cardio':'v','PRScr_geno':'8','PRScr_all':'8','all':'v'}
	
		
		#instantiate files to download
	allFilesTemp = glob.glob(f'{scoresPath}/*mixed.prs.csv')
	
	#cardioFile = [f for f in allFiles if 'Cardio' in f]
	
	#no covariate only data used in binning
	allFilesTemp = [f for f in allFilesTemp if 'covar' not in f]
	#remove the separate features from all model
	allFilesTemp = [f for f in allFilesTemp if 'FromAll' not in f]
	

	for holdout in [False,True]:

		if holdout:
			combinedPRSBinned = pd.DataFrame()
			outputPRSFile = f'{scoresPath}/combinedPRSGroups.holdout.csv'
			allFiles = [f for f in allFilesTemp if 'holdout' in f]
			
		else:
			prevalenceDf = pd.DataFrame()
			outputPRSFile = f'{scoresPath}/combinedPRSGroups.csv'
			outputPrevalenceFile = f'{scoresPath}/combinedPrevalencePRSGroups.csv'
			outputOddsRatioFile = f'{scoresPath}/combinedORPRSGroups.csv'
			allFiles = [f for f in allFilesTemp if 'holdout' not in f]
			
		combinedDf,prsList = combine_prs_files(allFiles)
		
		
		#scaled prs columns
		if holdout:
			try: #catch the exception if validation set hasn't been processed
				scaled_holdout = scale_holdout_data_manually(combinedDf, training_stats)
				combinedDf = scaled_holdout.copy()
			except NameError:
				print('ensure validation set has been combined and processed ..')
				pass
				
		else:
			#get the columns to scale
			columns_to_scale = [col for col in combinedDf.columns if 'prs' in col]
			scaled_df = scale_data(combinedDf[columns_to_scale])
			scaled_df.columns = [f'scaled_{col}' for col in columns_to_scale]
			combinedDf = combinedDf.merge(scaled_df,left_index=True,right_index=True,how='left')
			
			#bin the validation set
			for prs in prsList:
				prevalence, dfBinned = bin_prs(combinedDf,f'scaled_prs_{prs}',prs)
				prevalence['marker'] = marker_dict[prs]
				prevalence['color'] = color_dict[prs]
				prevalence['model'] = prs
				combinedDf = combinedDf.merge(dfBinned,left_index=True,right_index=True,how='left')
				prevalenceDf = pd.concat([prevalence,prevalenceDf],ignore_index=True)

				
			#get the training stats to scale prs in holdout
			training_stats = calculate_training_statistics(combinedDf)
			
			prevalenceDf.to_csv(outputPrevalenceFile,index=False)
			#calculate the raw PRS for each PRScr
		
			create_optimized_prevalence_plot(prevalenceDf,figPath,'validation')
			create_qq_plot_groups(combinedDf,figPath)
			
			percentileOR = calculate_odds_ratio_for_prs(combinedDf,'scaled_prs')
			percentileOR.to_csv(outputOddsRatioFile,index=False)
			
			combinedDf.drop(columns=['size','alpha'],inplace=True)
			
			
		combinedDf.reset_index().to_csv(outputPRSFile,index=False)
		
	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description="Combining PRS calculations ....")
	parser.add_argument("--pheno_path",help="path to pheno directory")
	
	
	args = parser.parse_args()
	
	# Prefer command-line input if provided; fallback to env var
	pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
	print(f"[PYTHON] Reading from: {pheno_path}")

#	pheno = 'celiacDisease'
#	pheno_path = f'/Users/kerimulterer/prsInteractive/results/{pheno}'
	if not pheno_path:
		raise ValueError("You must provide a data pheno path via --pheno_folder or set the PHENO_PATH environment variable.")
	main(pheno_path)

	