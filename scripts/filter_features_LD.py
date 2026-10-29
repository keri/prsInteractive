#!/usr/bin/env python3

import pandas as pd
import argparse
import os



def filter_epi_features(epiFeatures,ld2,rank_feature):
	'''Prune epi pairs that are LD-redundant with another pair.

	Two pairs (A1,A2) and (B1,B2) are considered redundant if either:

	  Case 1 — shared SNP: the pairs share exactly one SNP and their partner
	    (non-shared) SNPs are in LD.
	    e.g. (SNP1,SNP2) and (SNP1,SNP3): redundant if SNP2 <-> SNP3 in LD.

	  Case 2 — no shared SNP: both pairs have within-pair LD AND cross-pair LD
	    matches either the anti-diagonal pattern (A1<->B2 AND B1<->A2) or the
	    diagonal pattern (A1<->B1 AND A2<->B2).

	LD is checked bidirectionally so plink SNP_A/SNP_B ordering does not matter.
	The lower-ranked pair (by rank_feature) is pruned.
	'''

	def ld_between(snp1, snp2):
		return ld2[
			((ld2['SNP_A'] == snp1) & (ld2['SNP_B'] == snp2)) |
			((ld2['SNP_A'] == snp2) & (ld2['SNP_B'] == snp1))
		].shape[0] > 0

	ldfeaturesToFilter = []

	epiFeaturesToCheck = epiFeatures['feature'].tolist()
	for i in range(len(epiFeaturesToCheck)-1):
		feature = epiFeaturesToCheck[i]

		# skip pair-i if none of its SNPs appear anywhere in the LD file
		ldTemp = ld2[(ld2['SNP_A'].isin(feature.split(','))) | (ld2['SNP_B'].isin(feature.split(',')))]
		if ldTemp.shape[0] == 0:
			continue

		featureToFilter = [feature]
		A1 = feature.split(',')[0]
		A2 = feature.split(',')[1]

		for j in range(i+1,len(epiFeaturesToCheck)):
			feature2 = epiFeaturesToCheck[j]
			B1 = feature2.split(',')[0]
			B2 = feature2.split(',')[1]

			# Case 1: shared SNP — only partner SNPs need to be in LD
			if A1 == B1:
				redundant = ld_between(A2, B2)
			elif A1 == B2:
				redundant = ld_between(A2, B1)
			elif A2 == B1:
				redundant = ld_between(A1, B2)
			elif A2 == B2:
				redundant = ld_between(A1, B1)
			else:
				# Case 2: no shared SNP — within-pair LD on both sides plus
				# cross-pair LD in anti-diagonal or diagonal pattern
				pair_i_ld = ld_between(A1, A2)
				pair_j_ld = ld_between(B1, B2)
				anti_diag = ld_between(A1, B2) and ld_between(B1, A2)
				diag      = ld_between(A1, B1) and ld_between(A2, B2)
				redundant = pair_i_ld and pair_j_ld and (anti_diag or diag)

			if redundant:
				print('pair ', feature, ' is in LD with ', feature2)
				featureToFilter.append(feature2)

		epiFeaturesTemp = epiFeatures[epiFeatures['feature'].isin(featureToFilter)]
		try:
			epiFeaturesToPrune = epiFeaturesTemp.sort_values([rank_feature],ascending=False)['feature'].tolist()[1:]
		except KeyError:
			print(f'{rank_feature} is not found in dataset ....')
			print(epiFeatures.head())
			pass

		for f in epiFeaturesToPrune:
			ldfeaturesToFilter.append(f)

	return(list(set(ldfeaturesToFilter)))

def filter_main_in_ld(ld2,featuresMain,rank_feature):
	#match to feature 1
	
	#if both in main, then filter SNP with the lowest beta coefficient
	mainFeaturesInLD = []
	
	#filter data for those in which SNP_A is in main features
	mainLD1 = featuresMain[(featuresMain['feature'].isin(ld2['SNP_A'].tolist()))]
	

	#of those main features that is in LD, check to see if SNP_B is also in main 
	#if so capture the one that has the lowest beta coefficient to be pruned out
	#do this one feature at a time
	for ldSnp in mainLD1['feature']:
		#get the SNP_B in LD with SNP_A
		ldTemp = ld2[ld2['SNP_A'] == ldSnp]
		try:
			#find the main features that are in LD with ldSnp
			mainLD2 = featuresMain[(featuresMain['feature'].isin(ldTemp['SNP_B'].tolist()+[ldSnp]))]
			if mainLD2.shape[0] > 1:
				featuresToPrune = mainLD2.sort_values([rank_feature],ascending=False)['feature'].tolist()[1:]
				for f in featuresToPrune:
					mainFeaturesInLD.append(f)
		except KeyError:
			print(f'{rank_feature} not found in features dataset ...')
			print(featuresMain.head())
			pass
#	for f in featuresToPrune:
#		mainFeaturesInLD.append(f)

	return(mainFeaturesInLD)

#def main(phenoPath,features_filter_ld,pre_post_association):
def main(pheno_data,features_filter_ld):
	'''run the LD command separately which generates the plink.ld'''

	ld = pd.read_csv(f'{pheno_data}/finalModel.ld',sep='\s+')
	ld2 = ld[ld['R2'] > .6]
#	features = pd.read_csv(f'{phenoPath}/scores/importantFeaturesPostShap.csv')
	features = pd.read_csv(features_filter_ld)
	
	#this determines if the LD analysis is done before or after glm association modelling 
	#the 
#	if pre_post_association == 'post':
#	output_file = 'featureScoresReducedFinalModel.csv'
#	else:
#		output_file = 'importantFeaturesForAssociationAnalysis.csv'
	
	#refactored code includes a shap_zscore column not present in thesis results which means it will be the featureScoresReducedFinalModel created post glm modelling
	file_root = features_filter_ld.split('.')[0]
	output_file = f'{file_root}.filteredLD.csv'
	
	try:
		file_root = features_filter_ld.split('.')[0]
		output_file = f'{file_root}.filteredLD.csv'
		if 'iteration' in features.columns:
			features['shap_zscore'] = 2
		#get the main and epi weigths separate as the separate epi+main model performed best
		rank_feature = 'shap_zscore'
		featuresMain = features[features['data_type'] == 'main'][['feature','shap_zscore']]
		featuresEpi = features[(features['data_type'] == 'epi') & (features['feature'].str.contains(','))][['feature','shap_zscore']]
#		features.to_csv(f'{file_root}.preLD.csv',index=False)
	except KeyError:
		#the LD is done post glm modelling in early stages
		input_file_root = '/'.join(pheno_data.split('/')[:-1])
		input_file = f'{input_file_root}/featureScoresReducedFinalModel.csv'
		features = pd.read_csv(input_file)
		rank_feature = [col for col in features.columns if 'coef' in col][0]
		features2 = features[['model','feature',rank_feature]]
		featuresMain = features2[features2['model'] == 'main'][['feature',rank_feature]]
		featuresEpi = features2[(features2['model'] == 'epi') & (features2['feature'].str.contains(','))][['feature',rank_feature]]
	
	
	
	#############################################################################################
	#                              LD FEATURES TO PRUNE FOR MAIN SNPS                           #
	#############################################################################################

	
	mainFeaturesToPrune = filter_main_in_ld(ld2,featuresMain,rank_feature)
	
	#############################################################################################
	#                              LD FEATURES TO PRUNE FOR EPI PAIRS                           #
	#############################################################################################
	
	
	epiFeatureToPruneinLD = filter_epi_features(featuresEpi,ld2,rank_feature)
	
	featuresToPrune = mainFeaturesToPrune + epiFeatureToPruneinLD
	
	featuresFinal = features[~features['feature'].isin(featuresToPrune)]
	
	
	featuresFinal.to_csv(output_file,index=False)
	
	#if LD is done post glm modelling save the original dataset preLD
#	if pre_post_association == 'post':
#		features.to_csv(f'{phenoPath}/scores/featureScoresReducedFinalModel.preLD.csv',index=False)
	
	
		
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description="creating LD snp list ....")
	parser.add_argument("--pheno_data", help="Path to epi data folder (or pheno root for default LD run)")
	parser.add_argument("--features_to_ld_file", help="file used to filter ld")

	args = parser.parse_args()

	pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
	print(f"[PYTHON] Reading from: {pheno_data}")

	features_to_ld_file = args.features_to_ld_file or os.environ.get("FEATURES_TO_FILTER_LD")
	print(f'[PYTHON] Reading from feature files to filter for ld {features_to_ld_file}')

	if not pheno_data:
		raise ValueError("You must provide a data path via --pheno_data or set the PHENO_DATA environment variable.")

	if not features_to_ld_file:
		raise ValueError("You must provide a feature file to perform LD via --features_to_ld_file or set the FEATURES_TO_FILTER_LD environment variable.")

	main(pheno_data,features_to_ld_file)
	