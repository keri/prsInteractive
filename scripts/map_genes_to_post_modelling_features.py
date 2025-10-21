import pandas as pd
import numpy as np
import glob

def map_genes_to_features(featureDf,geneDf):
	'''merge genes and features into one dataframe'''
	########################  GENES MAPPED TO THE SNPS ###############
	matchingDf2 = featureDf.merge(geneDf[['feature','gene_names_feature1','gene_names_feature2','feature1','feature2']],on=['feature'],how='left')
	
	#Step 1: Sort the DataFrame based on the 'cohort' column
	matchingDf2 = matchingDf2.sort_values(by=['model','feature']).reset_index(drop=True)
	
	#fill feature 1 na when HLA loci
	matchingDf2.loc[matchingDf2['feature1'].isna(),'feature1'] = matchingDf2['feature']
	
	
	#fill na for snps that don't have a gene name
	matchingDf2.loc[matchingDf2['gene_names_feature1'].isna(),'gene_names_feature1'] = matchingDf2['feature1']
	
	
	#fill na for snp2 in pair that don't have a gene name
	matchingDf2.loc[((matchingDf2['gene_names_feature2'].isna()) & (matchingDf2['feature'].str.contains(','))),'gene_names_feature2'] = matchingDf2['feature2']
	
	
	#clean the -DT, -AS etc from gene names
	# Split the string only if it contains the substring, otherwise leave as is
	matchingDf2['gene_names_feature1'] = matchingDf2['gene_names_feature1'].apply(lambda x: x.split('-')[0] if pd.notna(x) and '-' in x else x)
	matchingDf2['gene_names_feature2'] = matchingDf2['gene_names_feature2'].apply(lambda x: x.split('-')[0] if pd.notna(x) and '-' in x else x)
	
#	
#	#shorten the ENSG to a truncated name
#	matchingDf2['gene_names_feature1'] = matchingDf2['gene_names_feature1'].apply(lambda x: 'E'+x[x.count('0')+4:] if pd.notna(x) and 'ENSG' in x else x)
#	matchingDf2['gene_names_feature2'] = matchingDf2['gene_names_feature2'].apply(lambda x: 'E'+x[x.count('0')+4:] if pd.notna(x) and 'ENSG' in x else x)
	
	
	
	# Combine columns with "," if the second column is not NaN
	matchingDf2['combined_genes'] = matchingDf2.apply(
		lambda row: f"{row['gene_names_feature1']},{row['gene_names_feature2']}" if pd.notna(row['gene_names_feature2']) else row['gene_names_feature1'],
		axis=1
	)
	return(matchingDf2)

def separate_cardio_features(df):
	#separate cardio features from geno features in cardio cohort
	df['cardioFeature'] = df.apply(lambda x: x['feature'].split(',')[0] if ((x['model'] == 'cardio') & ("," in x['feature'])) else '', axis=1)
	df['feature'] = df.apply(lambda x: ','.join(x['feature'].split(',')[1:]) if x['cardioFeature'] != '' else x['feature'], axis=1)
	
	return(df)


	
#def main(pheno):
pheno = 'type2Diabetes'
filePath = '/Users/kerimulterer/ukbiobank'
featurePath = f'{filePath}/{pheno}/tanigawaSet/featureScores/featureWeights'

featureFile = f"{featurePath}/featureScoresReducedFinalModelFilteredLD.csv"

#genes df with feature matched to genes
genes = pd.read_csv(f'{filePath}/{pheno}/tanigawaSet/networks/snpToProtein/combined{pheno}genes.csv')

if pheno == 'celiacDisease':
	knownAssociations = pd.read_csv(f'{filePath}/{pheno}/results/knownAssociations/gwas-association-downloaded_2025-02-16-EFO_0001060.tsv',sep="\t")
	knownAssociations = knownAssociations[['FIRST AUTHOR', 'DATE', 'JOURNAL','MAPPED_GENE','SNPS','OR or BETA', '95% CI (TEXT)']]
	#find the HLA region to merge separately due to different format
	dqKnownAssociations = knownAssociations[~knownAssociations['MAPPED_GENE'].isna() & knownAssociations['MAPPED_GENE'].str.contains('HLA')]
	snp_col = 'SNPS'
	knownAssociations.rename(columns={'OR or BETA':'previous_OR'},inplace=True)
else:
	snp_col = 'rsID'
	knownAssociations = pd.DataFrame()
	files = glob.glob(f'{filePath}/{pheno}/results/knownAssociations/*.txt')
	for file in files:
		df = pd.read_csv(file,sep='\t',usecols=['rsID','effect_weight','OR'])
		df = df[['rsID','effect_weight','OR']]
		df['JOURNAL'] = file.split('/')[-1].split('_')[0]
		df['95% CI (TEXT)'] = 'not_listed'
		knownAssociations = pd.concat([knownAssociations,df],ignore_index=True)
		
	#drop duplicates, keeping the OR rsID with the largest OR
	knownAssociations.sort_values(by=['rsID','OR'],ascending=False,inplace=True)
	knownAssociations.drop_duplicates(subset=['rsID'],keep='first',inplace=True)
	knownAssociations.rename(columns={'OR':'previous_OR'},inplace=True)



df = pd.read_csv(featureFile)
covariates = ['(Intercept)','age','SEX','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']

df = separate_cardio_features(df)
df = df[~df['feature'].isin(covariates)]
df = df[df['model'].isin(['cardio','main','epi','epi+main'])]


df2 = map_genes_to_features(df,genes)

#annotate loci data for a separate merge
df2['HLA_loci'] = df2['feature'].apply(lambda x: "HLA-" + x.split("_")[0] if 'D' in x else None)

#merge with gene df on first feature
merged1 = df2.merge(knownAssociations[~knownAssociations[snp_col].isna()],left_on='feature1',right_on=snp_col,how='left')

###################### get the SNP and HLA associations that match previously associated SNPs  ############

if pheno == 'celiacDisease':
	
	# Identify rows that didn't merge (NaN in df2 columns)
	unmatched = merged1[merged1[snp_col].isna()].drop(columns=knownAssociations.columns)

	# Merge unmatched rows on HLA_loci
	merged2 = unmatched.merge(knownAssociations[~knownAssociations['MAPPED_GENE'].isna()], left_on='HLA_loci', right_on='MAPPED_GENE', how='left')
	
	unmatched = merged2[merged2['MAPPED_GENE'].isna()].drop(columns=knownAssociations.columns)
	
	############### merge mapped gene for those SNPs that don't overlap with known SNPs  ########
	
	mergedGenes = unmatched.merge(knownAssociations[~knownAssociations['MAPPED_GENE'].isna()], left_on='gene_names_feature1', right_on='MAPPED_GENE', how='left')
	
	# Combine successfully merged rows from the first and second merge with the mapped genes merge
	df3 = pd.concat([merged1[~merged1['SNPS'].isna()], merged2[~merged2['MAPPED_GENE'].isna()], mergedGenes])
	df3.drop_duplicates(inplace=True)
	
else:
	df3 = merged1.copy()
	
df3['previously_known1'] = 'no'
df3.loc[~df3[snp_col].isna(),'previously_known1'] = 'yes'
	
####################  MERGE WITH SECOND FEATURE  ################

#merge with gene df on SECOND feature
merged2 = df2.merge(knownAssociations[~knownAssociations[snp_col].isna()],left_on='feature2',right_on=snp_col,how='left')

if pheno == 'celiacDisease':

	#merge on genes for those that do not match with SNPs
	# Identify rows that didn't merge (NaN in df2 columns)
	unmatched = merged2[merged2[snp_col].isna()].drop(columns=knownAssociations.columns)
	
	#there are no HLA genes in second feature so just need to merge with known genes
	mergedGenes = unmatched.merge(knownAssociations[~knownAssociations['MAPPED_GENE'].isna()], left_on='gene_names_feature2', right_on='MAPPED_GENE', how='left')
	
	df4 = pd.concat([merged2[~merged2[snp_col].isna()], mergedGenes])
	
else:
	df4 = merged2.copy()
	
#instantiate known associations column for second feature to combine in the end
df4['previously_known2'] = 'no'
df4.loc[~df4[snp_col].isna(),'previously_known2'] = 'yes'

df4.loc[df4['feature2'].isna(),'previously_known2'] = ' '
	
df4 = df4[['feature','model','lasso_coefs','feature2','gene_names_feature2','previously_known2'] + knownAssociations.columns.tolist()].merge(df3.drop(columns=['feature2','gene_names_feature2']),on=['feature','model','lasso_coefs'],suffixes=['_feature2','_feature1'],how='left')



df4['previously_known'] = df4['previously_known1'] + ','+ df4['previously_known2']
#df4['previous_OR'] = df4['OR_feature1'] + ',' + df4['OR_feature2']
#df4['previous_OR_CI'] = df4['95% CI (TEXT)_feature1'] + ','+ df4['95% CI (TEXT)_feature2']
#df4['journal'] = df4['JOURNAL_feature1'] + ','+ df4['JOURNAL_feature2']
df4.drop(columns=['previously_known1','previously_known2'],inplace=True)

df4.drop_duplicates(subset=['feature','model','lasso_coefs'],inplace=True)
df4['OR'] = np.exp(df4['lasso_coefs'])

#change model to G, GxG, G+GxG and GxGxE
model_dict = {'main':'G','epi':'GxG','epi+main':'G+(GxG)','cardio':'GxGxE'}
df4['PRS_calculation'] = df4['model'].apply( lambda x : model_dict[x])


df4 = df4[['feature', 'cardioFeature','PRS_calculation', 'lasso_coefs', 'OR','previously_known','feature1','gene_names_feature1','previous_OR_feature1', 'JOURNAL_feature1',
	'95% CI (TEXT)_feature1','feature2','gene_names_feature2', 'previous_OR_feature2','JOURNAL_feature2', '95% CI (TEXT)_feature2']]

#remove the trailing , from previously_known columns
df4['previously_known'] = df4['previously_known'].apply(lambda x : x.replace(', ',''))

df4 = df4[df4['lasso_coefs'] != 0]
df4.sort_values(['PRS_calculation','OR'],ascending=False,inplace=True)


df4.to_csv(f'{filePath}/{pheno}/results/supplementaryDataTables/combinedGeneFeaturePreviousOR.csv',index=False)