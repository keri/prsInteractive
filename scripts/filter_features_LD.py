#!/usr/bin/env python3

import pandas as pd
import argparse
import os



def filter_epi_features(epiFeatures,ld2,rank_feature):
	'''epiLD is pruned to the feature that is in LD based on plink output with all of ld and feature data from original LD match
	   epiLD columns = [feature,lasso_coefs_grid_search,feature1,feature2,CHR_A,BP_A,LD_SNP,CHR_B,BP_B,R2]
		features3 is the feature weights based on the final model containing both epi and main features
		features3 columns = [feature,lasso_coefs_grid_search,feature1,feature2]
		ldDf is plink output columns = [CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2]
		
		
		To determine if pair1(A1, B1) can tag pair2(A2, B2), perform the following calculations:
		
		Calculate LD Between Individual SNPs:
		
		Compute r² between A1 and A2
		
		Compute r² between A1 and B2
		
		Compute r² between B1 and A2
		
		Compute r² between B1 and B2
		
		Evaluate Pairwise LD:
		
		Assess whether the r² values meet your threshold (.6 for my purposes)
		
		Decision:
		
		If all computed r² values exceed the threshold, then (A1, B1) can be considered to tag (A2, B2) and only one needs to be kept.
		
		'''
	ldfeaturesToFilter = []
	
	
	epiFeaturesToCheck = epiFeatures['feature'].tolist()
	for i in range(len(epiFeaturesToCheck)-1):
		feature = epiFeaturesToCheck[i]
		
		#check to see if A1 or A2 are in LD with ANYTHING and if not can skip next step
		ldTemp = ld2[(ld2['SNP_A'].isin(feature.split(','))) | (ld2['SNP_B'].isin(feature.split(',')))]
		
		if ldTemp.shape[0] > 0: #if feature is in LD with something 
			#instantiate feature list with 1st feature
			featureToFilter = [feature]
			A1 = feature.split(',')[0]
			A2 = feature.split(',')[1]
			
			#compare against every other feature beginning at i+1
			for j in range(i+1,len(epiFeaturesToCheck)):
				ldNum = 0
				feature2 = epiFeaturesToCheck[j]
				B1 = feature2.split(',')[0]
				B2 = feature2.split(',')[1]
				#check to see if there is LD between A1 and A2
				#####################################   USE FOR STRICTER LD  #############################################
#				ldTemp2 = ldTemp[(ldTemp['SNP_A'].isin(feature2.split(','))) | (ldTemp['SNP_B'].isin(feature2.split(',')))]
				#check if A1 is in LD with A2
				ldTemp2 = ldTemp[((ldTemp['SNP_A'] == A1) & (ldTemp['SNP_B'] == A2))]
				if ldTemp2.shape[0] > 0:
					ldfeaturesToFilter.append(feature)
					ldNum += 1
				ldTemp2 = ldTemp[((ldTemp['SNP_A'] == A1) & (ldTemp['SNP_B'] == B2))]
				if ldTemp2.shape[0] > 0:
					ldNum += 1
				ldTemp2 = ldTemp[((ldTemp['SNP_A'] == B1) & (ldTemp['SNP_B'] == A2))]
				if ldTemp2.shape[0] > 0:
					ldNum += 1
				ldTemp2 = ldTemp[((ldTemp['SNP_A'] == B1) & (ldTemp['SNP_B'] == B2))]
				if ldTemp2.shape[0] > 0:
					ldNum += 1
#				if ldTemp2.shape[0] > 1:
				if ldNum == 4:
					print('pair ',feature, ' is in LD with ',feature2)
					featureToFilter.append(feature2)
					
			#get the tag pair with the highest coefficient
			epiFeaturesTemp = epiFeatures[epiFeatures['feature'].isin(featureToFilter)]
			try:
				epiFeaturesToPrune = epiFeaturesTemp.sort_values([rank_feature],ascending=False)['feature'].tolist()[1:]
			except KeyError:
				print(f'{rank_feature} is not found in dataset ....')
				print(epiFeatures.head())
				pass
				
			for f in epiFeaturesToPrune:
				ldfeaturesToFilter.append(f)
		else:
			pass
		
	
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
def main(phenoPath,features_filter_ld):
	'''run the LD command separately which generates the plink.ld'''
	
#	filePath = f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet'
	ld = pd.read_csv(f'{phenoPath}/finalModel.ld',sep='\s+')
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
		#get the main and epi weigths separate as the separate epi+main model performed best
		rank_feature = 'shap_zscore'
		featuresMain = features[features['data_type'] == 'main'][['feature','shap_zscore']]
		featuresEpi = features[(features['data_type'] == 'epi') & (features['feature'].str.contains(','))][['feature','shap_zscore']]
#		features.to_csv(f'{file_root}.preLD.csv',index=False)
	except KeyError:
		#the LD is done post glm modelling in early stages
		input_file_root = '/'.join(phenoPath.split('/')[:-1])
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
	parser.add_argument("--pheno_path", help="Path to the input pheno folder")
	parser.add_argument("--features_to_ld_file", help="file used to filter ld")
	
	
	
	args = parser.parse_args()
	
	# Prefer command-line input if provided; fallback to env var
#	pheno_path = '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
#	pre_post_association='pre'
	pheno_path = args.pheno_path or os.environ.get("PHENO_PATH")
	print(f"[PYTHON] Reading from: {pheno_path}")

#	pre_post_association = args.pre_post_association or os.environ.get("PRE_POST_ASSOCIATION")
#	print(f"pre or post association set to : {pre_post_association}")
	
	features_to_ld_file = args.features_to_ld_file or os.environ.get("FEATURES_TO_FILTER_LD")
	print(f'[PYTHON] Reading from feature files to filter for ld {features_to_ld_file}')
	
	#f"{pheno_data}/scores/importantFeaturesPostShap.csv"
#	print(f"[PYTHON] Reading features for LD from: {features_to_ld_file}")
	
#	pheno_path = '/Users/kerimulterer/prsInteractive/results/celiacDisease'
#	features_to_ld_file = f"{pheno_path}/productEpi/scores/importantFeaturesPostShap.csv"
	
	
	if not pheno_path:
		raise ValueError("You must provide a data pheno path via --pheno_path or set the PHENO_DATA environment variable.")
	
	if not features_to_ld_file:
		raise ValueError("You must provide a feature file to perform LD via --features_to_ld_file or set the FEATURES_TO_FILTER_LD environment variable.")	

		
	main(pheno_path,features_to_ld_file)
	