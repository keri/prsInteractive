#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
from scipy.stats import pearsonr
# sys.path.append("/Users/kerimulterer/ukbiobank/batchScripts/")
from helper.draw_plots import *
from helper.download import *
from helper.data_wrangling import *

def calculate_rsquared(df):
	corr, p_value = pearsonr(df['scaled_prs_main'], df['scaled_prs_epi'])
	return(corr,p_value)

def create_qq_plot_groups(combinedPRS,figurePath,use_all=use_all,use_epi_main=use_epi_main):
	''' use input to color cases in different colors based on group:
	input : dataframe columns = ['IID', 'PHENOTYPE', 'prs_epiMain', 'scaled_prs_epiMain', 'prs_main',
		'scaled_prs_main', 'prs_epi', 'scaled_prs_epi', 'use_epiMain',
		'use_epi', 'use_main', 'color', 'bin_epi+main', 'bin_main', 'bin_epi',
		'scaled_prs_final_group', 'scaled_ungrouped_optimized_prs',
		'ungrouped_bin_Final', 'grouped_optimized_prs',
		'scaled_grouped_optimized_prs', 'grouped_bin_Final']
	output : correlation plot with groups colors'''
	combinedPRS['size'] = 40
	combinedPRS['alpha'] = .6
	cases = combinedPRS[combinedPRS['PHENOTYPE'] == 2]
	cases['color'] = 'black'
	controls = combinedPRS[combinedPRS['PHENOTYPE'] == 1]
	controls['color'] = '#776cb1'
	

	######################## DEFINE GROUPS DYNAMICALLY ##################
	prsToUse = [f for f in combinedPRS.columns if 'scaled_prs' in f and 'prscr' not in f]
	
	if not use_epi_main:
		prsToUse = [f for f in prsToUse if 'epi+main' not in f]
		
	if not use_all:
		prsToUse = [f for f in prsToUse if 'all' in f]
		text_to_add = 'fromAll'
	
	# Build dfDict dynamically
	dfDict = {}
	use_cols = []
	for file in prsToUse:
		data_type = file.split('_')[-1]  # Extract just filename without path
		use_cols.append(f'bin_{data_type}')
		df = pd.read_csv(file)
		columns_to_drop = [col for col in df.columns if 'prs' in col]
		df.drop(columns=['PHENOTYPE']+columns_to_drop,inplace=True)
		dfDict[data_type] = df
	
	# Build conditions dynamically from whatever PRS types exist
	bin_columns = [f'bin_{data_type}' for data_type in dfDict.keys()]
	
	# Create masks by combining conditions
	high_risk_mask = pd.Series(False, index=groupedDf.index)
	low_risk_mask = pd.Series(False, index=groupedDf.index)
	
	for col in bin_columns:
		high_risk_mask |= (groupedDf[col] > 8)
		low_risk_mask |= (groupedDf[col] < 3)
		
	# Apply phenotype filters
	trainingDataHigh = groupedDf[high_risk_mask & (groupedDf['PHENOTYPE'] == 2)]
	trainingDataLow = groupedDf[low_risk_mask & (groupedDf['PHENOTYPE'] == 1)]

	str_text='combinedPRS'
	cases.loc[cases['bin_cardio'] > 8,'color'] = '#f2de19'
	cases.loc[cases['bin_epi+main'] > 8 ,'color'] = '#39c700'
	cases.loc[cases['bin_epi'] > 8,'color'] = '#19f0f2'
	cases.loc[cases['bin_main'] > 8,'color'] = '#C90016'
	cases.loc[cases['bin_all'] > 8,'color'] = '#c4771f'
		
	if use_all:
		str_text_final = f'{str_text}.WithAll'
		cases.loc[((cases['bin_main'] > 8) | (cases['bin_epi'] > 8) | (cases['bin_cardio'] > 8) | (cases['bin_epi+main'] > 8) | (cases['bin_all'] > 8)),'size'] = 60
		cases.loc[((cases['bin_main'] > 8) | (cases['bin_epi'] > 8) | (cases['bin_cardio'] > 8) | (cases['bin_epi+main'] > 8) | (cases['bin_all'] > 8)),'alpha'] = 1
	else:			
		str_text_final = str_text
		cases.loc[((cases['bin_main'] > 8) | (cases['bin_epi'] > 8) | (cases['bin_cardio'] > 8) | (cases['bin_epi+main'] > 8)),'size'] = 60
		cases.loc[((cases['bin_main'] > 8) | (cases['bin_epi'] > 8) | (cases['bin_cardio'] > 8) | (cases['bin_epi+main'] > 8)),'alpha'] = 1

	cases.loc[cases['color'] == '#f5a142','alpha'] = .6
	
	corr, p_value = pearsonr(combinedPRS['scaled_prs_main'], combinedPRS['scaled_prs_epi'])

	draw_plot(cases,controls,str_text,figurePath,corr,p_value)
			
			
def draw_validation_plot(df,figurePath,studies,validation_type='PGS'):
	#dfEpi is from PRSCombined with columns ['IID','scaled_prs_epi','bin_epi']
	#combinedPRS is the long form scaled PRS 
	
	combinedPRS = df.copy()
	
	#create long form 
#	combinedPRS = dfCopy.pivot(index=['IID','PHENOTYPE'], columns='PGS', values='scaled_prs').reset_index()
#	combinedPRS = dfCopy.merge(dfEpi,on=['IID'],how='left')
	combinedPRS['color'] = 'black'
	combinedPRS['size'] = 40
	combinedPRS['alpha'] = .5
	

	
	# calculate bins for each of the studies
	for g in studies:
		combinedPRS[f'bin_{g}'] = pd.qcut(combinedPRS[f'scaled_prs_{g}'], 10, labels=list(range(1,11)),duplicates='drop')
		#assign color to anyone with a bin > 8 in any study as one color
		combinedPRS.loc[(combinedPRS[f'bin_{g}'] > 8) ,'color'] = '#FFFF00'
		combinedPRS.loc[(combinedPRS[f'bin_{g}'] > 8) ,'size'] = 100
		combinedPRS.loc[(combinedPRS[f'bin_{g}'] > 8) ,'alpha'] = .8
		
	
	cases = combinedPRS[combinedPRS['PHENOTYPE'] == 2]
	controls = combinedPRS[combinedPRS['PHENOTYPE'] == 1]
	controls['color'] = '#776cb1'
	
	
	# find the cases that would be identified with at least ONE PRS calculation
	#color cases identified with other and not ours pop with yellow
	n_identified = 	cases.loc[cases['color'] == '#FFFF00'].shape[0]
	print('number of cases f that would be found with at least one PRS calculation', n_identified)
	
	#keep the color of people in main to the original color in our methods to compare against new methods
	str_text='allMainOtherVOurs.Cases'

	############   find the prs_main and prs_epi from our study that would be missed otherwise and found with all  ######
	#make the exclusive prs main pop with red
	cases.loc[((cases['bin_PGS000020'] < 9) & (cases['bin_PGS001357'] < 9) & (cases['bin_PGSVIC'] > 8)),'color'] = '#FF0000'
	cases.loc[cases['color'] == '#FF0000','size'] = 120
	cases.loc[cases['color'] == '#FF0000','alpha'] = 1
	n_missed = 	cases.loc[cases['color'] == '#FF0000'].shape[0]
#	cases.loc[cases['bin_PGSVIC'] > 8,'color'] = '#FF0000'
	print('number of main cases in our study missed with all other PRS ', n_missed)
	
	#color each study missed with all others
	#make the exclusive prs PGS000020 pop with green
	cases.loc[((cases['bin_PGS000020'] > 8) & (cases['bin_PGS001357'] < 9) & (cases['bin_PGSVIC'] < 9)),'color'] = '#008000'
	cases.loc[cases['color'] == '#008000','size'] = 100
	cases.loc[cases['color'] == '#008000','alpah'] = .8
	n_missed = 	cases.loc[cases['color'] == '#008000'].shape[0]
#	cases.loc[cases['bin_PGSVIC'] > 8,'color'] = '#FF0000'
	print('number of cases in 000020 study missed with all other PRS ', n_missed)
	
	#make the exclusive prs PGS001357 pop with yellow
	cases.loc[((cases['bin_PGS000020'] < 9) & (cases['bin_PGS001357'] > 8) & (cases['bin_PGSVIC'] < 9)),'color'] = '#f5a755'
	cases.loc[cases['color'] == '#f5a755','size'] = 100
	cases.loc[cases['color'] == '#f5a755','alpah'] = .8
	n_missed = 	cases.loc[cases['color'] == '#f5a755'].shape[0]
#	cases.loc[cases['bin_PGSVIC'] > 8,'color'] = '#FF0000'
	print('number of cases in 001357 study missed with all other PRS ', n_missed)
	
	
#	#make epi exclusive pop with blue
#	cases.loc[((cases['bin_PGS000020'] < 9) & (cases['bin_PGS001357'] < 9) & (cases['bin_PGSVIC'] < 9)),'color'] = '#0000FF'
#	cases.loc[cases['color'] == '#0000FF','size'] = 120
#	cases.loc[cases['color'] == '#0000FF','alpah'] = 1
#	n_missed = 	cases.loc[cases['color'] == '#0000FF'].shape[0]
#	print('number of epi cases in our study missed with all other PRS ', n_missed)
	
	#color cases identified using any and all PRS calculations pop with flourescent green
	cases.loc[((cases['bin_PGS000020'] > 8) & (cases['bin_PGS001357'] > 8) & (cases['bin_PGSVIC'] > 8)),'color'] = '#00FF00'
	cases.loc[cases['color'] == '#00FF00','size'] = 100
	cases.loc[cases['color'] == '#00FF00','alpha'] = 1
	n_missed = 	cases.loc[cases['color'] == '#00FF00'].shape[0]
	print('number of cases in our study that would be identified using all 3 PRS ', n_missed)
	
#	#color cases identified with other and not ours pop with yellow
#	cases.loc[(((cases['bin_PGS000020'] > 8) | (cases['bin_PGS001357'] > 8)) & ((cases['bin_PGSVIC'] < 9) & (cases['bin_epi'] < 9))),'color'] = '#fce303'
#	n_identified = 	cases.loc[cases['color'] == '#fce303'].shape[0]
#	print('number of cases f that would be found with other PRS and not identified with PRS (G) and PRS (GxG)', n_identified)
	
	#find number not identified with ANY PRS
	n_missed = cases.loc[((cases['bin_PGS000020'] < 9) & (cases['bin_PGS001357'] < 9) & (cases['bin_PGSVIC'] < 9))].shape[0]
	print('number of cases that would be not be identified with ANY PRS', n_missed)
	
	
	################ DRAW PLOT ########
	

	fig,ax = plt.subplots(figsize=(10,10))
		
	ax.scatter(y=controls['scaled_prs_PGSVIC'].values,x=controls['scaled_prs_PGS001357'].values,marker=".",c=controls['color'],s=controls['size'],alpha=controls['alpha'],label='Controls')
	ax.scatter(y=cases['scaled_prs_PGSVIC'].values,x=cases['scaled_prs_PGS001357'].values,marker="+",c=cases['color'],s=cases['size'],alpha=cases['alpha'],label='Cases')
		
	ax.set_xlabel('scaled PRS : PGS001357')
	ax.set_ylabel('scaled PRS : G')
	ax.grid(True) 
	ax.set_title(f'Standardized PRS across 3 studies : Main effect SNPs')
	handles,labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels,loc='upper left',fontsize=15)
#   rsquared = round((linreg.rvalue*linreg.rvalue),2)
	
	
	rsquared1,pvalue1 = pearsonr(combinedPRS['scaled_prs_PGSVIC'], combinedPRS['scaled_prs_PGS001357'])
	print('r squared prs main v PGS001357 :',rsquared1)
	print('pvalue prs main v PGS001357 :',pvalue1)
	
	
	rsquared2,pvalue2 = pearsonr(combinedPRS['scaled_prs_PGSVIC'], combinedPRS['scaled_prs_PGS000020'])
	print('r squared prs main v PGS000020 :',rsquared2)
	print('pvalue prs main v PGS000020 :',pvalue2)
	
	# Add text anchored to the bottom-right of the plot
#	plt.text(
#		1.0, .01,               # Position (relative to the axes)
#		f'pearson r-sqared : {rsquared}',        # Text to display
#		transform=plt.gca().transAxes,  # Use Axes coordinates (0, 0 is bottom-left, 1, 1 is top-right)
#		fontsize=12,             # Text size
#		verticalalignment='bottom',    # Align text vertically to the bottom
#		horizontalalignment='right'    # Align text horizontally to the right
#	)
#	
#	plt.text(
#		.7, .05,               # Position (relative to the axes)
#		f'p value : {pvalue}',        # Text to display
#		transform=plt.gca().transAxes,  # Use Axes coordinates (0, 0 is bottom-left, 1, 1 is top-right)
#		fontsize=12,             # Text size
#		verticalalignment='bottom',    # Align text vertically to the bottom
#		horizontalalignment='right'    # Align text horizontally to the right
#	)
	
	plt.savefig(f'{figurePath}/{str_text}.coloredStudies.QQColorPlotMainOnly.png')
	return(combinedPRS)
	
	
			

def draw_plot(cases,controls,str_text,figurePath,rsquared,pvalue):
	fig,ax = plt.subplots(figsize=(10,10))
#	model_type = f'{str_text}.{suffix}'
#   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
	ax.scatter(x=controls['scaled_prs_epi'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'],s=controls['size'],alpha=controls['alpha'],label='Controls')
#   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
	ax.scatter(x=cases['scaled_prs_epi'].values,y=cases['scaled_prs_main'].values,marker="+",c=cases['color'],s=cases['size'],alpha=cases['alpha'],label='Cases')
#   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
	
	ax.set_xlabel('epi only')
	ax.set_ylabel('main only')
	ax.grid(True) 
	ax.set_title(f'Standardized {str_text} : Main V Epi')
	handles,labels = ax.get_legend_handles_labels()
	ax.legend(handles,labels,loc='upper left',fontsize=15)
#   rsquared = round((linreg.rvalue*linreg.rvalue),2)
	
	# rsquared,pvalue = calculate_rsquared(pd.concat([cases,controls]))
	# Add text anchored to the bottom-right of the plot
	plt.text(
		1.0, .01,               # Position (relative to the axes)
		f'pearson r-sqared : {rsquared}',        # Text to display
		transform=plt.gca().transAxes,  # Use Axes coordinates (0, 0 is bottom-left, 1, 1 is top-right)
		fontsize=12,             # Text size
		verticalalignment='bottom',    # Align text vertically to the bottom
		horizontalalignment='right'    # Align text horizontally to the right
	)
	
	plt.text(
		.7, .05,               # Position (relative to the axes)
		f'p value : {pvalue}',        # Text to display
		transform=plt.gca().transAxes,  # Use Axes coordinates (0, 0 is bottom-left, 1, 1 is top-right)
		fontsize=12,             # Text size
		verticalalignment='bottom',    # Align text vertically to the bottom
		horizontalalignment='right'    # Align text horizontally to the right
	)
	
	plt.savefig(f'{figurePath}/{str_text}.QQColorPlot.png')
	ax.clear()
	cases['color'] = 'black'
	
def main(phenoData,use_all=use_all,use_epi_main=use_epi_main):
	scoresPath = f'{phenoData}/scores'
	figPath = f'{phenoData}/figures'

	combinedPRS = pd.read_csv(f'{scoresPath}/combinedPRSGroups.csv')
#	if pheno == 'type2Diabetes':
	create_qq_plot_groups(combinedPRS,figPath,use_epi_main=use_epi_main)
		

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description="calculating model performance of G vs other using predictions of trained models and PRS calculations....")
	parser.add_argument("--pheno_path", help="Path to the input pheno folder")
	parser.add_argumtn("--use_all",default=False,help="set as true if QQ plot will include PRS calculated from modelled combined features (default=False)")
	parser.add_argumtn("--use_epi_main",default=False,help="set as true if QQ plot will include PRS calculated from combined epi+main models (default=False)")
	
	args = parser.parse_args()
	
	# Prefer command-line input if provided; fallback to env var
	pheno_data = args.pheno_path or os.environ.get("PHENO_DATA")
	print(f"[PYTHON] Reading from: {pheno_data}")
	
	use_all = args.use_all or os.environ.get("USE_ALL")
	print(f"[PYTHON] setting use_all to: {use_all}")
	
	use_epi_main = args.use_epi_main or os.environ.get("USE_EPI_MAIN")
	print(f"[PYTHON] setting use_epi_main to: {use_epi_main}")
	
	#pheno_data = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi'
	
	
	if not pheno_data:
		raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
		
	
	main(pheno_data,use_all=use_all,use_epi_main=use_epi_main)