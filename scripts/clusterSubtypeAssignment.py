#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
from matplotlib.colors import LinearSegmentedColormap


# ============================================================================
# CALCULATE FUNCTIONAL CLUSTER RANKING ENRICHMENT
# ============================================================================

SUBTYPE_COLORS = {
	'SAID': '#e28743',
	'SIDD': '#873e23',
	'SIRD': '#eeeee4',
	'MOD': '#21130d',
	'MARD': '#eab676',
	'UNKNOWN': '#938369'
}

COHORT_COLORS = {
	'GxGxE': '#f2de19',
	'GxG': '#19f0f2',
	'G+(GxG)': '#39c700',
	'G': '#C90016'
}

def expand_cohort_rows(df,cohort_col='cohort'):
	"""
	Provides expansion of cohort rows for inclusion of genes important to more than one OR exclusively
	Parameters:
	-----------
	df : pandas.DataFrame
		DataFrame with gene/feature data
	cohort_col : str
		Name of the column containing cohort information (can be comma-separated)
		
	Returns:
	--------
	df with:
		- 'expanded_df': DataFrame of filtered/expanded rows (cohorts × clusters)

	"""
	
	# Identify functional cluster columns
	cluster_cols = [col for col in df.columns if col.endswith('_count')]
	cluster_cols = [col for col in cluster_cols if 'opathy' not in col]
	cluster_cols = [col for col in cluster_cols if 'wnt' not in col]
	cluster_cols = [col for col in cluster_cols if 'mitochondria' not in col]
	
	
	# Expand rows where cohort column contains multiple comma-separated cohorts
	expanded_rows = []
	for _, row in df.iterrows():
		cohort_value = row[cohort_col]
		if pd.isna(cohort_value):
			continue
		
		# Split by comma and strip whitespace
		cohorts = [c.strip() for c in str(cohort_value).split(',')]
		
		#for cohort in cohorts:
		# Create a row for each cohort
		if len(cohorts) == 1:
			#if cohort:  # Skip empty strings
			row_copy = row.copy()
			row_copy[cohort_col] = cohorts[0]
			expanded_rows.append(row_copy)
			
	# Create expanded dataframe
	df_expanded = pd.DataFrame(expanded_rows)
	
	return df_expanded,cluster_cols
	

def rank_functional_clusters_by_cohort(df, cohort_col='cohort', 
											or_col='OR', 
											use_effect_size=True):
	"""
	Provides normalized ranking of functional clusters per cohort.
	Handles multiple cohorts per gene (comma-separated values).
	
	Parameters:
	-----------
	df : pandas.DataFrame
		DataFrame with gene/feature data
	cohort_col : str
		Name of the column containing cohort information (can be comma-separated)
	or_col : str
		Name of the column containing odds ratios
	weight_by_or : bool
		If True, weight cluster associations by odds ratio
		
	Returns:
	--------
	dict with:
		- 'scores_df': DataFrame of normalized scores (cohorts × clusters)
		- 'cluster_ranks_df': DataFrame of cluster ranks within each cohort
		- 'cohort_ranks_df': DataFrame of cohort ranks for each cluster
		- 'top_clusters': DataFrame showing top cluster per cohort
		- 'top_cohorts': DataFrame showing top cohort per cluster
		- 'summary': Dict with analysis metadata
	"""
	
#	# Identify functional cluster columns
#	cluster_cols = [col for col in df.columns if col.endswith('_count')]
#	cluster_cols = [col for col in cluster_cols if 'opathy' not in col]
#	cluster_cols = [col for col in cluster_cols if 'wnt' not in col]
#	cluster_cols = [col for col in cluster_cols if 'mitochondria' not in col]
#	
#	
#	# Expand rows where cohort column contains multiple comma-separated cohorts
#	expanded_rows = []
#	for _, row in df.iterrows():
#		cohort_value = row[cohort_col]
#		if pd.isna(cohort_value):
#			continue
#		
#		# Split by comma and strip whitespace
#		cohorts = [c.strip() for c in str(cohort_value).split(',')]
#		
#		#for cohort in cohorts:
#		# Create a row for each cohort
#		if len(cohorts) == 1:
#			#if cohort:  # Skip empty strings
#			row_copy = row.copy()
#			row_copy[cohort_col] = cohorts[0]
#			expanded_rows.append(row_copy)
#				
#	# Create expanded dataframe
#	df_expanded = pd.DataFrame(expanded_rows)
	df_expanded,cluster_cols = expand_cohort_rows(df, 'cohort')
		
	# Initialize results
	cohort_scores = {}
	cohort_gene_counts = {}
	cohort_directions = {}  # Track risk vs protective
	cohort_protective_scores = {}
	cohort_risk_scores = {}
	
	# Calculate scores for each cohort
	for cohort in df_expanded[cohort_col].unique():
		cohort_data = df_expanded[df_expanded[cohort_col] == cohort]
		n_genes = len(cohort_data)
		cohort_gene_counts[cohort] = n_genes
		
		cluster_scores = {}
		cluster_directions = {}
		protective_scores = {}
		risk_scores = {}
		
		for cluster in cluster_cols:
			total_score = 0
			total_protective = 0
			total_risk = 0
			direction_counts = {'risk': 0, 'protective': 0, 'neutral': 0}
			
			for _, row in cohort_data.iterrows():
				# Check if associated with this cluster
				if pd.notna(row.get(cluster)) and row[cluster] == 1:
					or_value = row.get(or_col, 1)
					
					if use_effect_size:
						# Use effect size (treats risk and protective equally)
						effect_size, direction = calculate_or_effect_size(or_value)
						total_score += effect_size
						
						# Track separately
						if direction == 'protective':
							total_protective += effect_size
						elif direction == 'risk':
							total_risk += effect_size
					else:
						# Original behavior (raw OR)
						total_score += or_value if pd.notna(or_value) else 1
						
						if or_value < 1:
							direction = 'protective'
							total_protective += or_value
						elif or_value > 1:
							direction = 'risk'
							total_risk += or_value
						else:
							direction = 'neutral'
							
					direction_counts[direction] += 1
					
			# Normalize by number of genes
			cluster_scores[cluster] = total_score / n_genes if n_genes > 0 else 0
			protective_scores[cluster] = total_protective / n_genes if n_genes > 0 else 0
			risk_scores[cluster] = total_risk / n_genes if n_genes > 0 else 0
			
			# Determine dominant direction
			if direction_counts['risk'] > direction_counts['protective']:
				cluster_directions[cluster] = 'risk'
			elif direction_counts['protective'] > direction_counts['risk']:
				cluster_directions[cluster] = 'protective'
			else:
				cluster_directions[cluster] = 'mixed'
				
		cohort_scores[cohort] = cluster_scores
		cohort_directions[cohort] = cluster_directions
		cohort_protective_scores[cohort] = protective_scores
		cohort_risk_scores[cohort] = risk_scores
		
	# Convert to DataFrames
	scores_df = pd.DataFrame(cohort_scores).T
	scores_df.index.name = 'cohort'
	
	directions_df = pd.DataFrame(cohort_directions).T
	directions_df.index.name = 'cohort'
	
	protective_df = pd.DataFrame(cohort_protective_scores).T
	protective_df.index.name = 'cohort'
	
	risk_df = pd.DataFrame(cohort_risk_scores).T
	risk_df.index.name = 'cohort'
	
	# Rank clusters within each cohort (1 = strongest effect)
	cluster_ranks_df = scores_df.rank(axis=1, ascending=False, method='min')
	
	# Rank cohorts for each cluster (1 = strongest effect)
	cohort_ranks_df = scores_df.rank(axis=0, ascending=False, method='min')
	
	# Identify top cluster per cohort
	top_cluster_names = scores_df.idxmax(axis=1)
	top_cluster_scores = scores_df.max(axis=1)

	# Get directions for top clusters (pandas 2.0 compatible)
	top_cluster_directions = [
		directions_df.loc[cohort, cluster] 
		for cohort, cluster in zip(scores_df.index, top_cluster_names)
	]

	top_clusters = pd.DataFrame({
		'top_cluster': top_cluster_names.str.replace('_count', ''),
		'score': top_cluster_scores,
		'direction': top_cluster_directions,
		'n_genes': pd.Series(cohort_gene_counts)
	})

	# Identify top cohort per cluster
	top_cohorts = pd.DataFrame({
		'top_cohort': scores_df.idxmax(axis=0),
		'score': scores_df.max(axis=0)
	})
	top_cohorts.index = top_cohorts.index.str.replace('_count', '')
	top_cohorts.index.name = 'cluster'

	# Summary statistics
	summary = {
		'total_cohorts': len(cohort_scores),
		'total_clusters': len(cluster_cols),
		'total_genes_original': len(df),
		'total_genes_expanded': len(df_expanded),
		'genes_added_by_splitting': len(df_expanded) - len(df),
		'genes_per_cohort': cohort_gene_counts,
		'use_effect_size': use_effect_size
	}

	return {
		'scores_df': scores_df,
		'cluster_ranks_df': cluster_ranks_df,
		'cohort_ranks_df': cohort_ranks_df,
		'top_clusters': top_clusters,
		'top_cohorts': top_cohorts,
		'effect_directions': directions_df,
		'protective_scores': protective_df,
		'risk_scores': risk_df,
		'summary': summary
	}
		
	
def display_direction_analysis(results):
	"""
	Display analysis of risk vs protective effects.
	
	Parameters:
	-----------
	results : dict
		Results from rank_functional_clusters_by_cohort_corrected
	"""
	
	print("\n" + "=" * 80)
	print("RISK vs PROTECTIVE EFFECTS ANALYSIS")
	print("=" * 80)
	
	directions_df = results['effect_directions']
	
	# Count direction types per cohort
	print("\nCluster directions by cohort:")
	print("-" * 80)
	
	for cohort in directions_df.index:
		direction_counts = directions_df.loc[cohort].value_counts()
		print(f"\n{cohort}:")
		for direction, count in direction_counts.items():
			print(f"  {direction}: {count} clusters")
			
	# Show top protective clusters
	print("\n" + "=" * 80)
	print("TOP PROTECTIVE CLUSTERS (OR < 1)")
	print("=" * 80)
	
	protective_df = results['protective_scores']
	
	for cohort in protective_df.index:
		top_protective = protective_df.loc[cohort].nlargest(3)
		top_protective = top_protective[top_protective > 0]
		
		if len(top_protective) > 0:
			print(f"\n{cohort}:")
			for cluster, score in top_protective.items():
				cluster_name = cluster.replace('_count', '')
				direction = results['effect_directions'].loc[cohort, cluster]
				if direction in ['protective', 'mixed']:
					print(f"  {cluster_name}: {score:.4f} (protective effect)")
					
	# Show top risk clusters
	print("\n" + "=" * 80)
	print("TOP RISK CLUSTERS (OR > 1)")
	print("=" * 80)
	
	risk_df = results['risk_scores']
	
	for cohort in risk_df.index:
		top_risk = risk_df.loc[cohort].nlargest(3)
		top_risk = top_risk[top_risk > 0]
		
		if len(top_risk) > 0:
			print(f"\n{cohort}:")
			for cluster, score in top_risk.items():
				cluster_name = cluster.replace('_count', '')
				direction = results['effect_directions'].loc[cohort, cluster]
				if direction in ['risk', 'mixed']:
					print(f"  {cluster_name}: {score:.4f} (risk effect)")
					

#def display_results(results, top_n=3):
#	"""Display analysis results in a readable format"""
#	
#	print("=" * 80)
#	print("SUMMARY STATISTICS")
#	print("=" * 80)
#	print(f"Total cohorts: {results['summary']['total_cohorts']}")
#	print(f"Total clusters: {results['summary']['total_clusters']}")
#	print(f"Total genes/features: {results['summary']['total_genes']}")
#	print(f"Weighted by OR: {results['summary']['weight_by_or']}")
#	
#	print("\n" + "=" * 80)
#	print("TOP CLUSTER PER COHORT (normalized by number of genes)")
#	print("=" * 80)
#	print(results['top_clusters'].sort_values('score', ascending=False).to_string())
#	
#	print("\n" + "=" * 80)
#	print(f"TOP {top_n} COHORTS PER CLUSTER")
#	print("=" * 80)
#	
#	# For each cluster, show top cohorts
#	scores_df = results['scores_df']
#	for cluster in scores_df.columns[:5]:  # Show first 5 clusters as example
#		cluster_name = cluster.replace('_count', '')
#		print(f"\n{cluster_name}:")
#		top_cohorts = scores_df[cluster].nlargest(top_n)
#		for cohort, score in top_cohorts.items():
#			n_genes = results['summary']['genes_per_cohort'][cohort]
#			print(f"  {cohort} (n={n_genes}): {score:.4f}")
#			
#	print("\n" + "=" * 80)
#	print("NORMALIZED SCORES (sample - showing first 5 cohorts × 5 clusters)")
#	print("=" * 80)
#	print(results['scores_df'].iloc[:5, :5].round(4).to_string())
#	
#	print("\n" + "=" * 80)
#	print("CLUSTER RANKINGS PER COHORT (sample - 1 = top cluster)")
#	print("=" * 80)
#	print(results['cluster_ranks_df'].iloc[:5, :5].to_string())

def calculate_weighted_subtype_enrichment(df, results, cluster_subtype_mapping,
																					cohort_col='cohort',
																					or_col='OR',
																					weighting='effect_size'):
	"""
	Calculate T2D subtype enrichment weighted by effect size and optionally distinctiveness.
	
	This determines what proportion of each subtype's EFFECT (not just gene count)
	comes from each cohort, accounting for OR magnitude.
	
	Parameters:
	-----------
	df : pandas.DataFrame
			Original data with genes and cluster associations
	results : dict
			Results from rank_functional_clusters_by_cohort_corrected
	cluster_subtype_mapping : DataFrame
			Mapping from load_cluster_subtype_mapping()
	cohort_col : str
			Name of cohort column
	or_col : str
			Name of OR column
	weighting : str
			'effect_size': Weight by log(OR) magnitude
			'distinctiveness': Weight by effect size × distinctiveness
			'normalized': Weight by normalized OR contribution
			
	Returns:
	--------
	dict with:
			- 'cohort_weights_per_subtype': Weighted contributions
			- 'cohort_percentages_per_subtype': Percentage of subtype effect from each cohort
			- 'gene_counts': Raw gene counts for comparison
			- 'direction_breakdown': Risk vs protective breakdown
	"""
	
		# Expand rows with comma-separated cohorts
#		expanded_rows = []
#		for _, row in df.iterrows():
#				cohort_value = row[cohort_col]
#				if pd.isna(cohort_value):
#						continue
#			
#				cohorts = [c.strip() for c in str(cohort_value).split(',')]
#			
#				for cohort in cohorts:
#						if cohort:
#								row_copy = row.copy()
#								row_copy[cohort_col] = cohort
#								expanded_rows.append(row_copy)
#							
#		df_expanded = pd.DataFrame(expanded_rows)
#	
#		# Get cluster columns
#		cluster_cols = [col for col in df_expanded.columns if col.endswith('_count')]
	
	df_expanded,cluster_cols = expand_cohort_rows(df, 'cohort')

	# Create mapping dictionary
	cluster_to_subtypes = {}
	for cluster in cluster_cols:
			cluster_name = cluster.replace('_count', '')
			matching = cluster_subtype_mapping[
					cluster_subtype_mapping['functional_cluster'] == cluster_name
			]
			if len(matching) > 0:
					cluster_to_subtypes[cluster] = matching['subtype'].tolist()
			else:
					cluster_to_subtypes[cluster] = []
				
	# Get distinctiveness scores if using that weighting
	distinctiveness = None
	if weighting == 'distinctiveness' and 'scores_df' in results:
			scores_df = results['scores_df']
			z_scores = scores_df.apply(lambda col: (col - col.mean()) / col.std(), axis=0)
			percentile_ranks = scores_df.rank(axis=0, pct=True) * 100
		
			distinctiveness = pd.DataFrame(index=scores_df.index, columns=scores_df.columns)
			for cluster in scores_df.columns:
					for cohort in scores_df.index:
							score = scores_df.loc[cohort, cluster]
							z_score = z_scores.loc[cohort, cluster]
							percentile = percentile_ranks.loc[cohort, cluster]
						
							if score > 0:
									distinctiveness.loc[cohort, cluster] = abs(z_score * percentile / 100)
							else:
									distinctiveness.loc[cohort, cluster] = 0
								
			distinctiveness = distinctiveness.astype(float)
		
	# Calculate weighted contributions
	subtype_cohort_weights = {}
	subtype_cohort_counts = {}
	subtype_cohort_risk = {}
	subtype_cohort_protective = {}

	for cohort in df_expanded[cohort_col].unique():
			cohort_data = df_expanded[df_expanded[cohort_col] == cohort]
		
			for _, row in cohort_data.iterrows():
					for cluster in cluster_cols:
							if pd.notna(row.get(cluster)) and row[cluster] == 1:
									subtypes = cluster_to_subtypes.get(cluster, [])
									or_value = row.get(or_col, 1)
								
									# Calculate effect size
									if pd.notna(or_value) and or_value > 0:
											effect_size = abs(np.log(or_value))
											direction = 'risk' if or_value > 1 else 'protective' if or_value < 1 else 'neutral'
									else:
											effect_size = 0
											direction = 'neutral'
										
									# Apply weighting
									if weighting == 'effect_size':
											weight = effect_size
									elif weighting == 'distinctiveness':
											# Combine effect size with distinctiveness
											if distinctiveness is not None:
													dist_score = distinctiveness.loc[cohort, cluster]
													weight = effect_size * (1 + dist_score)  # Boost by distinctiveness
											else:
													weight = effect_size
									elif weighting == 'normalized':
											# Use normalized score from results
											if 'scores_df' in results:
													weight = results['scores_df'].loc[cohort, cluster]
											else:
													weight = effect_size
									else:
											weight = 1  # Equal weighting
										
									# Count contributions
									for subtype in subtypes:
											# Initialize dictionaries
											if subtype not in subtype_cohort_weights:
													subtype_cohort_weights[subtype] = {}
													subtype_cohort_counts[subtype] = {}
													subtype_cohort_risk[subtype] = {}
													subtype_cohort_protective[subtype] = {}
												
											if cohort not in subtype_cohort_weights[subtype]:
													subtype_cohort_weights[subtype][cohort] = 0
													subtype_cohort_counts[subtype][cohort] = 0
													subtype_cohort_risk[subtype][cohort] = 0
													subtype_cohort_protective[subtype][cohort] = 0
												
											# Add weighted contribution
											subtype_cohort_weights[subtype][cohort] += weight
											subtype_cohort_counts[subtype][cohort] += 1
										
											# Track direction
											if direction == 'risk':
													subtype_cohort_risk[subtype][cohort] += weight
											elif direction == 'protective':
													subtype_cohort_protective[subtype][cohort] += weight
												
	# Convert to DataFrames
	weights_df = pd.DataFrame(subtype_cohort_weights).T.fillna(0)
	weights_df.index.name = 'subtype'

	counts_df = pd.DataFrame(subtype_cohort_counts).T.fillna(0).astype(int)
	counts_df.index.name = 'subtype'

	risk_df = pd.DataFrame(subtype_cohort_risk).T.fillna(0)
	risk_df.index.name = 'subtype'

	protective_df = pd.DataFrame(subtype_cohort_protective).T.fillna(0)
	protective_df.index.name = 'subtype'

	# Calculate percentages
	percentages_df = weights_df.div(weights_df.sum(axis=1), axis=0) * 100

	# Summary
	enrichment_summary = {
			'total_subtypes': len(weights_df),
			'cohorts_contributing': weights_df.columns.tolist(),
			'weighting_method': weighting,
			'total_weighted_contribution_per_subtype': weights_df.sum(axis=1).to_dict()
	}

	return {
			'cohort_weights_per_subtype': weights_df,
			'cohort_percentages_per_subtype': percentages_df,
			'gene_counts': counts_df,
			'risk_contributions': risk_df,
			'protective_contributions': protective_df,
			'enrichment_summary': enrichment_summary
	}
		

def plot_weighted_subtype_pies(subtype_results, save_path='weighted_subtype_composition.png'):
	"""
	Create pie charts showing weighted cohort composition within each subtype.
	
	Parameters:
	-----------
	subtype_results : dict
			Results from calculate_weighted_subtype_enrichment
	save_path : str
			Path to save the figure
	"""
	

	percentages_df = subtype_results['cohort_percentages_per_subtype']
	weights_df = subtype_results['cohort_weights_per_subtype']
	counts_df = subtype_results['gene_counts']

	# Custom colors
	SUBTYPE_COLORS = {
			'SAID': '#e28743',
			'SIDD': '#873e23',
			'SIRD': '#eeeee4',
			'MOD': '#21130d',
			'MARD': '#eab676',
			'UNKNOWN': '#938369'
	}

	COHORT_COLORS = {
			'GxGxE': '#f2de19',
			'GxG': '#19f0f2',
			'G+(GxG)': '#39c700',
			'G': '#C90016'
	}

	subtypes = percentages_df.index.tolist()
	n_subtypes = len(subtypes)

	# Create figure
	if n_subtypes <= 5:
			fig, axes = plt.subplots(1, n_subtypes, figsize=(5 * n_subtypes, 6))
	else:
			n_cols = 3
			n_rows = int(np.ceil(n_subtypes / n_cols))
			fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 6 * n_rows))
		
	if n_subtypes == 1:
			axes = np.array([axes])
	axes = axes.flatten()

	# Create pie chart for each subtype
	for idx, subtype in enumerate(subtypes):
		ax = axes[idx]
	
		# Get cohort contributions
		data = percentages_df.loc[subtype]
		data = data[data > 0].sort_values(ascending=False)
	
		if len(data) > 0:
			# Get colors
			pie_colors = [COHORT_COLORS.get(c, '#95a5a6') for c in data.index]
		
			# Create pie chart
#			wedges, texts, autotexts = ax.pie(
#				data.values,
#				labels=data.index,
#				colors=pie_colors,
#				autopct='%1.1f%%',
#				startangle=90,
#				pctdistance=0.85
#			)
		
			wedges, texts = ax.pie(
				data.values,
				labels=data.index,
				colors=pie_colors,
				startangle=90,
				pctdistance=0.85
			)
		
			# Improve text visibility
			for text in texts:
				text.set_fontsize(9)
				text.set_fontweight('bold')
#				for autotext in autotexts:
#					autotext.set_color('white')
#					autotext.set_fontweight('bold')
#					autotext.set_fontsize(8)
				
			#Add title with weighted contribution and gene count
			total_weight = weights_df.loc[subtype].sum()
			total_genes = counts_df.loc[subtype].sum()
			title_color = SUBTYPE_COLORS.get(subtype, '#2d3748')
		
			ax.set_title(
				f'{subtype}', 
				fontsize=13, fontweight='bold', color=title_color, pad=20
			)
		
			# Add note about weighting
			method = subtype_results['enrichment_summary']['weighting_method']
			ax.text(0, -1.3, f'Weighted by: {method}', 
							ha='center', fontsize=8, style='italic', color='gray')
		else:
			ax.text(0.5, 0.5, f'{subtype}\nNo data', 
							ha='center', va='center', fontsize=12)
			ax.axis('off')
				
	# Hide unused subplots
	for idx in range(n_subtypes, len(axes)):
			axes[idx].axis('off')
		
	# Add overall title
	fig.suptitle('Cohort Contribution to Each T2D Subtype\n(Weighted by Effect Size)', 
							fontsize=16, fontweight='bold', y=0.98)

	plt.tight_layout()
	plt.savefig(save_path, dpi=300, bbox_inches='tight')
	print(f"\n✓ Weighted subtype pie charts saved to: {save_path}")

	return fig


def plot_contribution_comparison(subtype_results, save_path='contribution_comparison.png'):
	"""
	Compare raw counts vs weighted contributions for each subtype.
	
	Parameters:
	-----------
	subtype_results : dict
			Results from calculate_weighted_subtype_enrichment
	save_path : str
			Path to save the figure
	"""

	counts_df = subtype_results['gene_counts']
	weights_df = subtype_results['cohort_weights_per_subtype']

	# Get dominant cohort by each method
	dominant_by_count = counts_df.idxmax(axis=1)
	dominant_by_weight = weights_df.idxmax(axis=1)

	# Create comparison
	comparison = pd.DataFrame({
			'Subtype': counts_df.index,
			'By_Count': dominant_by_count.values,
			'By_Weight': dominant_by_weight.values,
			'Count_%': [counts_df.loc[st, dominant_by_count[st]] / counts_df.loc[st].sum() * 100 
									for st in counts_df.index],
			'Weight_%': [weights_df.loc[st, dominant_by_weight[st]] / weights_df.loc[st].sum() * 100 
										for st in weights_df.index]
	})

	# Identify where they differ
	comparison['Different'] = comparison['By_Count'] != comparison['By_Weight']

	# Create figure
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

	# Plot 1: Bar comparison
	x = np.arange(len(comparison))
	width = 0.35

	bars1 = ax1.bar(x - width/2, comparison['Count_%'], width, 
									label='By Gene Count', color='#3498db', alpha=0.8)
	bars2 = ax1.bar(x + width/2, comparison['Weight_%'], width, 
									label='By Effect Weight', color='#e74c3c', alpha=0.8)

	ax1.set_xlabel('Subtype', fontsize=12, fontweight='bold')
	ax1.set_ylabel('Dominant Cohort Contribution (%)', fontsize=12, fontweight='bold')
	ax1.set_title('Gene Count vs Effect-Weighted Contribution', fontsize=14, fontweight='bold')
	ax1.set_xticks(x)
	ax1.set_xticklabels(comparison['Subtype'], rotation=45, ha='right')
	ax1.legend()
	ax1.grid(True, alpha=0.3, axis='y')

	# Highlight differences
	for i, diff in enumerate(comparison['Different']):
			if diff:
					ax1.text(i, max(comparison.loc[i, 'Count_%'], comparison.loc[i, 'Weight_%']) + 2,
									'⚠', ha='center', fontsize=16, color='orange')
				
	# Plot 2: Show which cohort is dominant
	methods = ['By Count', 'By Weight']
	subtypes = comparison['Subtype'].tolist()

	# Create grid showing dominant cohort
	dominant_matrix = np.zeros((len(subtypes), 2))
	cohort_names = sorted(set(comparison['By_Count'].tolist() + comparison['By_Weight'].tolist()))
	cohort_to_idx = {c: i for i, c in enumerate(cohort_names)}

	for i, row in comparison.iterrows():
			dominant_matrix[i, 0] = cohort_to_idx[row['By_Count']]
			dominant_matrix[i, 1] = cohort_to_idx[row['By_Weight']]
		
	im = ax2.imshow(dominant_matrix, cmap='tab10', aspect='auto')

	ax2.set_xticks([0, 1])
	ax2.set_xticklabels(methods)
	ax2.set_yticks(np.arange(len(subtypes)))
	ax2.set_yticklabels(subtypes)
	ax2.set_title('Dominant Cohort by Method', fontsize=14, fontweight='bold')

	# Add text annotations
	for i in range(len(subtypes)):
			for j in range(2):
					cohort_idx = int(dominant_matrix[i, j])
					cohort = cohort_names[cohort_idx]
					text = ax2.text(j, i, cohort, ha="center", va="center", 
												color="white", fontweight='bold', fontsize=9)
				
	# Highlight differences with border
	for i, diff in enumerate(comparison['Different']):
			if diff:
					rect = plt.Rectangle((0-0.5, i-0.5), 2, 1, fill=False, 
															edgecolor='orange', linewidth=3)
					ax2.add_patch(rect)
				
	plt.tight_layout()
	plt.savefig(save_path, dpi=300, bbox_inches='tight')
	print(f"✓ Contribution comparison saved to: {save_path}")

	return fig, comparison


def display_weighted_results(subtype_results):
	"""
	Display weighted subtype enrichment results.
	
	Parameters:
	-----------
	subtype_results : dict
			Results from calculate_weighted_subtype_enrichment
	"""

	print("\n" + "=" * 80)
	print("WEIGHTED SUBTYPE ENRICHMENT ANALYSIS")
	print("=" * 80)

	method = subtype_results['enrichment_summary']['weighting_method']
	print(f"\nWeighting method: {method}")
	print(f"Subtypes identified: {', '.join(subtype_results['cohort_weights_per_subtype'].index.tolist())}")

	print("\n" + "=" * 80)
	print("WEIGHTED CONTRIBUTIONS (Effect-Size Based)")
	print("=" * 80)
	print(subtype_results['cohort_weights_per_subtype'].round(3).to_string())

	print("\n" + "=" * 80)
	print("PERCENTAGE OF SUBTYPE EFFECT FROM EACH COHORT")
	print("=" * 80)
	print(subtype_results['cohort_percentages_per_subtype'].round(1).to_string())

	print("\n" + "=" * 80)
	print("RAW GENE COUNTS (For Comparison)")
	print("=" * 80)
	print(subtype_results['gene_counts'].to_string())

	# Compare dominant cohorts
	print("\n" + "=" * 80)
	print("DOMINANT COHORT COMPARISON")
	print("=" * 80)

	weights_df = subtype_results['cohort_weights_per_subtype']
	counts_df = subtype_results['gene_counts']
	percentages_df = subtype_results['cohort_percentages_per_subtype']

	for subtype in weights_df.index:
			dominant_weight = weights_df.loc[subtype].idxmax()
			dominant_count = counts_df.loc[subtype].idxmax()
		
			weight_pct = percentages_df.loc[subtype, dominant_weight]
		
			print(f"\n{subtype}:")
			print(f"  By effect weight: {dominant_weight} ({weight_pct:.1f}%)")
			print(f"  By gene count: {dominant_count} ({counts_df.loc[subtype, dominant_count]} genes)")
		
			if dominant_weight != dominant_count:
					print(f"  ⚠ Different dominant cohorts! Effect-weighted shows stronger contribution from {dominant_weight}")
				
				
	

def calculate_distinctive_clusters(results):
	"""
	Calculate which clusters are distinctively important to each cohort
	by comparing their scores to other cohorts.
	
	A cluster is "distinctive" if it scores high in one cohort but low in others.
	
	Returns:
	--------
	dict with:
		- 'z_scores': Z-score normalized scores (how many SDs above/below mean)
		- 'percentile_ranks': Percentile rank of each score across cohorts
		- 'distinctiveness': Combined distinctiveness score
		- 'top_distinctive': Top distinctive clusters per cohort
	"""
	scores_df = results['scores_df']
	
	# Method 1: Z-score normalization (for each cluster across cohorts)
	z_scores = scores_df.apply(lambda col: (col - col.mean()) / col.std(), axis=0)
	
	# Method 2: Percentile ranking (for each cluster across cohorts)
	percentile_ranks = scores_df.rank(axis=0, pct=True) * 100
	
	# Method 3: Distinctiveness score
	# High distinctiveness = high absolute score + high relative rank
	distinctiveness = pd.DataFrame(index=scores_df.index, columns=scores_df.columns)
	
	for cluster in scores_df.columns:
		for cohort in scores_df.index:
			score = scores_df.loc[cohort, cluster]
			z_score = z_scores.loc[cohort, cluster]
			percentile = percentile_ranks.loc[cohort, cluster]
			
			# Distinctiveness: combine z-score and percentile
			# High z-score means it's unusually high for this cluster
			# We weight by the actual score to ensure it's meaningfully high
			if score > 0:
				distinctiveness.loc[cohort, cluster] = z_score * percentile / 100
			else:
				distinctiveness.loc[cohort, cluster] = 0
				
	# Convert to numeric
	distinctiveness = distinctiveness.astype(float)
	
	# Find top distinctive clusters per cohort
	top_distinctive = pd.DataFrame(index=scores_df.index, columns=[
		'top_distinctive_cluster', 
		'distinctiveness_score',
		'z_score',
		'percentile',
		'absolute_score',
		'rank_in_cohort'
	])
	
	for cohort in scores_df.index:
		# Get cluster with highest distinctiveness
		top_cluster = distinctiveness.loc[cohort].idxmax()
		
		top_distinctive.loc[cohort] = {
			'top_distinctive_cluster': top_cluster.replace('_count', ''),
			'distinctiveness_score': distinctiveness.loc[cohort, top_cluster],
			'z_score': z_scores.loc[cohort, top_cluster],
			'percentile': percentile_ranks.loc[cohort, top_cluster],
			'absolute_score': scores_df.loc[cohort, top_cluster],
			'rank_in_cohort': results['cluster_ranks_df'].loc[cohort, top_cluster]
		}
		
	return {
		'z_scores': z_scores,
		'percentile_ranks': percentile_ranks,
		'distinctiveness': distinctiveness,
		'top_distinctive': top_distinctive
	}
	

def compare_cluster_across_cohorts(results, cluster_name):
	"""
	Compare a specific cluster's performance across all cohorts.
	
	Parameters:
	-----------
	results : dict
		Results from rank_functional_clusters_by_cohort
	cluster_name : str
		Name of cluster to analyze (with or without '_count' suffix)
	
	Returns:
	--------
	DataFrame with cluster scores, ranks, and statistics across cohorts
	"""
	if not cluster_name.endswith('_count'):
		cluster_name = cluster_name + '_count'
		
	scores_df = results['scores_df']
	cluster_ranks_df = results['cluster_ranks_df']
	
	if cluster_name not in scores_df.columns:
		raise ValueError(f"Cluster '{cluster_name}' not found")
		
	comparison = pd.DataFrame({
		'cohort': scores_df.index,
		'score': scores_df[cluster_name].values,
		'rank': cluster_ranks_df[cluster_name].values,
		'n_genes': [results['summary']['genes_per_cohort'][c] for c in scores_df.index]
	})
	
	# Add percentile and z-score
	comparison['percentile'] = (comparison['score'].rank(pct=True) * 100).round(1)
	comparison['z_score'] = ((comparison['score'] - comparison['score'].mean()) / 
							comparison['score'].std()).round(3)
	
	# Sort by score descending
	comparison = comparison.sort_values('score', ascending=False).reset_index(drop=True)
	
	return comparison


def find_cohort_signature_clusters(results, cohort_name, top_n=5):
	"""
	Find the most distinctive/signature clusters for a specific cohort.
	
	These are clusters that:
	1. Score relatively high in this cohort
	2. Score relatively low in other cohorts (making them distinctive)
	
	Parameters:
	-----------
	results : dict
		Results from rank_functional_clusters_by_cohort
	cohort_name : str
		Name of cohort to analyze
	top_n : int
		Number of top distinctive clusters to return
	
	Returns:
	--------
	DataFrame with signature clusters and their statistics
	"""
	if cohort_name not in results['scores_df'].index:
		raise ValueError(f"Cohort '{cohort_name}' not found")
		
	distinctive_results = calculate_distinctive_clusters(results)
	
	# Get distinctiveness scores for this cohort
	cohort_distinctiveness = distinctive_results['distinctiveness'].loc[cohort_name]
	cohort_scores = results['scores_df'].loc[cohort_name]
	cohort_ranks = results['cluster_ranks_df'].loc[cohort_name]
	cohort_z_scores = distinctive_results['z_scores'].loc[cohort_name]
	cohort_percentiles = distinctive_results['percentile_ranks'].loc[cohort_name]
	
	# Combine into DataFrame
	signature = pd.DataFrame({
		'cluster': cohort_distinctiveness.index,
		'distinctiveness': cohort_distinctiveness.values,
		'absolute_score': cohort_scores.values,
		'rank_in_cohort': cohort_ranks.values,
		'z_score': cohort_z_scores.values,
		'percentile_vs_other_cohorts': cohort_percentiles.values
	})
	
	# Clean cluster names
	signature['cluster'] = signature['cluster'].str.replace('_count', '')
	
	# Sort by distinctiveness and get top N
	signature = signature.sort_values('distinctiveness', ascending=False).head(top_n)
	signature = signature.reset_index(drop=True)
	
	return signature

def calculate_or_effect_size(or_value):
	"""
	Convert OR to effect size that treats risk and protective effects equally.
	
	Uses log transformation to make the scale symmetric:
	- OR = 2.0 (100% increase) → log(2) = 0.693
	- OR = 0.5 (50% decrease) → log(0.5) = -0.693 → abs = 0.693
	
	Parameters:
	-----------
	or_value : float
		Odds ratio value
	
	Returns:
	--------
	tuple: (effect_size, direction)
		- effect_size: Absolute magnitude of effect (always positive)
		- direction: 'risk', 'protective', or 'neutral'
	"""
	if pd.isna(or_value) or or_value <= 0:
		return 0, 'neutral'
	
	if or_value > 1:
		direction = 'risk'
	elif or_value < 1:
		direction = 'protective'
	else:
		direction = 'neutral'
		
	# Log transformation makes scale symmetric
	effect_size = abs(np.log(or_value))
	
	return effect_size, direction

	



# ============================================================================
# LOAD CLUSTER-SUBTYPE MAPPING
# ============================================================================

def load_cluster_subtype_mapping(mapping_file='clusterToSubtype.csv'):
	"""
	Load and expand cluster-to-subtype mapping.
	Handles clusters that map to multiple subtypes (comma-separated).
	
	Parameters:
	-----------
	mapping_file : str
		Path to CSV file with columns: functional_cluster, subtype
	
	Returns:
	--------
	DataFrame with columns:
		- functional_cluster: cluster name (e.g., 'immunity')
		- subtype: T2D subtype (e.g., 'SAID')
		- cluster_with_suffix: cluster name with '_count' (e.g., 'immunity_count')
	
	Example:
	--------
	>>> mapping = load_cluster_subtype_mapping('clusterToSubtype.csv')
	>>> print(mapping.head())
	"""
	# Read the mapping file
	mapping_df = pd.read_csv(mapping_file)
	
	#remove the UNKNOWN SUBTYPE
	mapping_df = mapping_df[mapping_df['functional_cluster'] != 'UNKNOWN']
	
	# Expand rows where subtype contains multiple comma-separated values
	expanded_mappings = []
	
	for _, row in mapping_df.iterrows():
		cluster = row['functional_cluster']
		subtypes_str = str(row['subtype'])
		
		# Split by comma and strip whitespace
		subtypes = [s.strip() for s in subtypes_str.split(',')]
		
		# Create a row for each subtype
		for subtype in subtypes:
			if subtype and subtype != 'nan':  # Skip empty or NaN values
				expanded_mappings.append({
					'functional_cluster': cluster,
					'subtype': subtype
				})
			
	# Convert to DataFrame
	expanded_df = pd.DataFrame(expanded_mappings)
	
	# Add '_count' suffix to match cluster column names in data
	expanded_df['cluster_with_suffix'] = expanded_df['functional_cluster'] + '_count'
	
	print(f"Loaded {len(mapping_df)} original mappings")
	print(f"Expanded to {len(expanded_df)} mappings (handling multi-subtype entries)")
	print(f"Subtypes found: {expanded_df['subtype'].unique().tolist()}")
	
	return expanded_df


# ============================================================================
# CALCULATE SUBTYPE ENRICHMENT
# ============================================================================

def calculate_subtype_enrichment(df, results, cluster_subtype_mapping, 
									cohort_col='cohort'):
	"""
	Calculate T2D subtype enrichment showing cohort composition within each subtype.
	
	This determines what proportion of each subtype's associations come from each cohort.
	
	Parameters:
	-----------
	df : pandas.DataFrame
		Original data with genes and cluster associations
	results : dict
		Results from rank_functional_clusters_by_cohort
	cluster_subtype_mapping : DataFrame
		Mapping from load_cluster_subtype_mapping()
	cohort_col : str
		Name of cohort column
		
	Returns:
	--------
	dict with:
		- 'cohort_counts_per_subtype': DataFrame (subtypes × cohorts) with counts
		- 'cohort_percentages_per_subtype': DataFrame (subtypes × cohorts) with percentages
		- 'enrichment_summary': Summary statistics
	
	Example:
	--------
	>>> mapping = load_cluster_subtype_mapping('clusterToSubtype.csv')
	>>> subtype_results = calculate_subtype_enrichment(df, results, mapping)
	>>> print(subtype_results['cohort_percentages_per_subtype'])
	"""
	
	# Expand rows with comma-separated cohorts (same as in main function)
#	expanded_rows = []
#	for _, row in df.iterrows():
#		cohort_value = row[cohort_col]
#		if pd.isna(cohort_value):
#			continue
#		
#		cohorts = [c.strip() for c in str(cohort_value).split(',')]
#		
#		for cohort in cohorts:
#			if cohort:
#				row_copy = row.copy()
#				row_copy[cohort_col] = cohort
#				expanded_rows.append(row_copy)
#				
#	df_expanded = pd.DataFrame(expanded_rows)
#	
#	# Get cluster columns
#	cluster_cols = [col for col in df_expanded.columns if col.endswith('_count')]
	
	df_expanded,cluster_cols = expand_cohort_rows(df, 'cohort')
	
	# Create mapping dictionary for quick lookup
	cluster_to_subtypes = {}
	for cluster in cluster_cols:
		cluster_name = cluster.replace('_count', '')
		matching = cluster_subtype_mapping[
			cluster_subtype_mapping['functional_cluster'] == cluster_name
		]
		if len(matching) > 0:
			cluster_to_subtypes[cluster] = matching['subtype'].tolist()
		else:
			cluster_to_subtypes[cluster] = []
			
	# Calculate subtype-cohort counts (what cohorts contribute to each subtype)
	subtype_cohort_counts = {}
	
	for cohort in df_expanded[cohort_col].unique():
		cohort_data = df_expanded[df_expanded[cohort_col] == cohort]
		
		# For each gene in the cohort
		for _, row in cohort_data.iterrows():
			# Check each cluster
			for cluster in cluster_cols:
				# If gene is associated with this cluster
				if pd.notna(row.get(cluster)) and row[cluster] == 1:
					# Get subtypes for this cluster
					subtypes = cluster_to_subtypes.get(cluster, [])
					
					# Count this cohort's contribution to each subtype
					for subtype in subtypes:
						if subtype not in subtype_cohort_counts:
							subtype_cohort_counts[subtype] = {}
						if cohort not in subtype_cohort_counts[subtype]:
							subtype_cohort_counts[subtype][cohort] = 0
						subtype_cohort_counts[subtype][cohort] += 1
						
	# Convert to DataFrame (rows = subtypes, columns = cohorts)
	counts_df = pd.DataFrame(subtype_cohort_counts).T.fillna(0).astype(int)
	counts_df.index.name = 'subtype'
	
	# Calculate percentages (what % of each subtype comes from each cohort)
	percentages_df = counts_df.div(counts_df.sum(axis=1), axis=0) * 100
	
	# Summary statistics
	enrichment_summary = {
		'total_subtypes': len(counts_df),
		'cohorts_contributing': counts_df.columns.tolist(),
		'total_associations_per_subtype': counts_df.sum(axis=1).to_dict()
	}
	
	return {
		'cohort_counts_per_subtype': counts_df,
		'cohort_percentages_per_subtype': percentages_df,
		'enrichment_summary': enrichment_summary
	}
	

	

# ============================================================================
# CREATE PIE CHARTS
# ============================================================================



def plot_signature_distinctiveness_scatter(results, save_path='signature_distinctiveness_scatter.png'):
	"""
	Create a scatter plot showing distinctiveness of top signature clusters.
	
	X-axis: Top signature clusters (ranking 1-2) for each cohort
	Y-axis: Distinctiveness value for all cohorts
	
	Points are color-coded by cohort, showing which clusters are truly
	distinctive (high y-value) and which cohorts they belong to.
	
	Parameters:
	-----------
	results : dict
			Results from rank_functional_clusters_by_cohort
	save_path : str
			Path to save the figure
	
	Returns:
	--------
	matplotlib figure object
	
	Example:
	--------
	>>> fig = plot_signature_distinctiveness_scatter(results)
	>>> plt.show()
	"""

	# Calculate distinctiveness
	scores_df = results['scores_df']

	# Z-score normalization for each cluster across cohorts
	z_scores = scores_df.apply(lambda col: (col - col.mean()) / col.std(), axis=0)

	# Percentile ranking
	percentile_ranks = scores_df.rank(axis=0, pct=True) * 100

	# Distinctiveness score
	distinctiveness = pd.DataFrame(index=scores_df.index, columns=scores_df.columns)
	for cluster in scores_df.columns:
			for cohort in scores_df.index:
					score = scores_df.loc[cohort, cluster]
					z_score = z_scores.loc[cohort, cluster]
					percentile = percentile_ranks.loc[cohort, cluster]
				
					if score > 0:
							distinctiveness.loc[cohort, cluster] = z_score * percentile / 100
					else:
							distinctiveness.loc[cohort, cluster] = 0
						
	distinctiveness = distinctiveness.astype(float)

	# Get cluster rankings within each cohort
	cluster_ranks = scores_df.rank(axis=1, ascending=False, method='min').astype(int)

	# Collect top 1-2 signature clusters for each cohort
	plot_data = []

	for cohort in scores_df.index:
			# Get clusters ranked 1 or 2
			top_clusters = cluster_ranks.loc[cohort][cluster_ranks.loc[cohort] <= 2].index.tolist()
		
			for cluster in top_clusters:
					rank = cluster_ranks.loc[cohort, cluster]
					dist_value = distinctiveness.loc[cohort, cluster]
					cluster_name = cluster.replace('_count', '')
				
					plot_data.append({
							'cohort': cohort,
							'cluster': cluster_name,
							'rank': rank,
							'distinctiveness': dist_value,
							'z_score': z_scores.loc[cohort, cluster],
							'absolute_score': scores_df.loc[cohort, cluster]
					})
				
	plot_df = pd.DataFrame(plot_data)

	# Create figure
	fig, ax = plt.subplots(figsize=(14, 8))

	# Plot each cohort with its own color
	cohorts = plot_df['cohort'].unique()

	for cohort in cohorts:
			cohort_data = plot_df[plot_df['cohort'] == cohort]
			color = COHORT_COLORS.get(cohort, '#95a5a6')
		
			# Separate rank 1 and rank 2 for different markers
			rank1 = cohort_data[cohort_data['rank'] == 1]
			rank2 = cohort_data[cohort_data['rank'] == 2]
		
			if len(rank1) > 0:
					ax.scatter(
							rank1['cluster'],
							rank1['distinctiveness'],
							s=200,
							c=[color],
							marker='o',
							alpha=0.8,
							edgecolors='black',
							linewidth=2,
							label=f'{cohort} (Rank 1)'
					)
				
			if len(rank2) > 0:
					ax.scatter(
							rank2['cluster'],
							rank2['distinctiveness'],
							s=150,
							c=[color],
							marker='s',
							alpha=0.6,
							edgecolors='black',
							linewidth=1.5,
							label=f'{cohort} (Rank 2)'
					)
				
	# Add horizontal line at distinctiveness = 0
	ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

	# Formatting
	ax.set_xlabel('Signature Cluster', fontsize=12, fontweight='bold')
	ax.set_ylabel('Distinctiveness Score', fontsize=12, fontweight='bold')
	ax.set_title('Top Signature Clusters by Distinctiveness\n(Ranks 1-2 per Cohort)', 
							fontsize=14, fontweight='bold', pad=20)

	# Rotate x-axis labels
	plt.xticks(rotation=45, ha='right')

	# Add legend
	ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, 
						fancybox=True, shadow=True, fontsize=9)

	# Add grid
	ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

	# Add annotation for interpretation
	textstr = ('◯ = Rank 1 (Top cluster)\n'
						'□ = Rank 2 (Second cluster)\n'
						'Higher Y = More distinctive')
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
					verticalalignment='top', bbox=props)

	plt.tight_layout()
	plt.savefig(save_path, dpi=300, bbox_inches='tight')
	print(f"✓ Signature distinctiveness scatter plot saved to: {save_path}")

	# Print summary
	print("\n" + "=" * 80)
	print("SIGNATURE DISTINCTIVENESS SUMMARY")
	print("=" * 80)

	# Find most distinctive signature clusters
	top_distinctive = plot_df.nlargest(5, 'distinctiveness')
	print("\nTop 5 most distinctive signature clusters:")
	for idx, row in top_distinctive.iterrows():
			print(f"  {row['cohort']} → {row['cluster']} (Rank {row['rank']}): "
						f"Distinctiveness = {row['distinctiveness']:.3f}")
		
	return fig
	
# ============================================================================
# PLOT RISK/PROTECTION RESULTS
# ============================================================================
	
def plot_risk_protective_comparison(results, save_path='risk_protective_comparison.png'):
	"""
	Create a comparison plot showing risk vs protective effects per cohort.
	
	Parameters:
	-----------
	results : dict
		Results from rank_functional_clusters_by_cohort_corrected
	save_path : str
		Path to save the figure
	"""
	
	protective_df = results['protective_scores']
	risk_df = results['risk_scores']
	
	# Sum protective and risk scores per cohort
	cohort_protective = protective_df.sum(axis=1).sort_values(ascending=False)
	cohort_risk = risk_df.sum(axis=1).sort_values(ascending=False)
	
	# Create comparison DataFrame
	comparison = pd.DataFrame({
		'Protective': cohort_protective,
		'Risk': cohort_risk
	})
	
	# Create figure
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
	
	# Plot 1: Stacked bar chart
	comparison.plot(kind='bar', stacked=False, ax=ax1, color=['#2ecc71', '#e74c3c'])
	ax1.set_title('Risk vs Protective Effects by Cohort', fontsize=14, fontweight='bold')
	ax1.set_xlabel('Cohort', fontsize=12)
	ax1.set_ylabel('Total Effect Size', fontsize=12)
	ax1.legend(title='Effect Type')
	ax1.grid(True, alpha=0.3, axis='y')
	plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
	
	# Plot 2: Risk-Protective Ratio
	comparison['Ratio'] = comparison['Risk'] / (comparison['Protective'] + 0.001)  # Avoid division by zero
	comparison['Ratio'].plot(kind='bar', ax=ax2, color='#3498db')
	ax2.axhline(y=1, color='gray', linestyle='--', linewidth=1, label='Equal balance')
	ax2.set_title('Risk/Protective Ratio by Cohort', fontsize=14, fontweight='bold')
	ax2.set_xlabel('Cohort', fontsize=12)
	ax2.set_ylabel('Risk/Protective Ratio', fontsize=12)
	ax2.legend()
	ax2.grid(True, alpha=0.3, axis='y')
	plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
	
	plt.tight_layout()
	plt.savefig(save_path, dpi=300, bbox_inches='tight')
	print(f"\n✓ Risk vs protective comparison saved to: {save_path}")
	
	return fig


# ============================================================================
# DISPLAY SUBTYPE RESULTS
# ============================================================================

def display_subtype_results(subtype_results):
	"""
	Display subtype enrichment results in readable format.
	
	Parameters:
	-----------
	subtype_results : dict
		Results from calculate_subtype_enrichment
	
	Example:
	--------
	>>> display_subtype_results(subtype_results)
	"""
	
	print("\n" + "=" * 80)
	print("T2D SUBTYPE ENRICHMENT ANALYSIS")
	print("=" * 80)
	
	print(f"\nSubtypes identified: {', '.join(subtype_results['cohort_counts_per_subtype'].index.tolist())}")
	print(f"Total cohorts analyzed: {len(subtype_results['cohort_counts_per_subtype'].columns)}")
	
	print("\n" + "=" * 80)
	print("COHORT COUNTS WITHIN EACH SUBTYPE")
	print("(How many associations from each cohort contribute to each subtype)")
	print("=" * 80)
	print(subtype_results['cohort_counts_per_subtype'].to_string())
	
	print("\n" + "=" * 80)
	print("COHORT PERCENTAGES WITHIN EACH SUBTYPE")
	print("(What % of each subtype comes from each cohort)")
	print("=" * 80)
	print(subtype_results['cohort_percentages_per_subtype'].round(1).to_string())
	
	# Find dominant cohort per subtype
	print("\n" + "=" * 80)
	print("DOMINANT COHORT FOR EACH SUBTYPE")
	print("=" * 80)
	
	for subtype, row in subtype_results['cohort_percentages_per_subtype'].iterrows():
		if row.sum() > 0:
			dominant = row.idxmax()
			percentage = row.max()
			total = subtype_results['cohort_counts_per_subtype'].loc[subtype].sum()
			print(f"{subtype}: {dominant} contributes {percentage:.1f}% (n={int(total)} total associations)")
			
	# Show subtype composition details
	print("\n" + "=" * 80)
	print("SUBTYPE COMPOSITION DETAILS")
	print("=" * 80)
	
	for subtype in subtype_results['cohort_percentages_per_subtype'].index:
		total = subtype_results['cohort_counts_per_subtype'].loc[subtype].sum()
		print(f"\n{subtype} (n={int(total)} associations):")
		
		composition = subtype_results['cohort_percentages_per_subtype'].loc[subtype].sort_values(ascending=False)
		for cohort, pct in composition.items():
			if pct > 0:
				count = subtype_results['cohort_counts_per_subtype'].loc[subtype, cohort]
				print(f"  {cohort}: {pct:.1f}% (n={int(count)})")

# ============================================================================
# CREATE HEATMAP
# ============================================================================
				
def parse_cluster_subtype_mapping(mapping_df: pd.DataFrame) -> Dict[str, List[str]]:
	"""
	Parse cluster-to-subtype mapping from a DataFrame.
	
	Expected DataFrame format (one row per cluster-subtype pair):
	functional_cluster | subtype | cluster_with_suffix
	beta_cell_expression | SIDD | beta_cell_expression_count
	obesity | MOD | obesity_count
	obesity | SIRD | obesity_count
	
	Parameters:
	-----------
	mapping_df : pd.DataFrame
		DataFrame with 'functional_cluster', 'subtype', and 'cluster_with_suffix' columns
	
	Returns:
	--------
	tuple: (cluster_to_subtypes_dict, cluster_to_column_dict)
		- cluster_to_subtypes_dict: Maps functional_cluster to list of subtypes
		- cluster_to_column_dict: Maps functional_cluster to cluster_with_suffix column name
	"""
	cluster_to_subtypes = {}
	cluster_to_column = {}
	
	for _, row in mapping_df.iterrows():
		cluster = row['functional_cluster']
		subtype = row['subtype']
		cluster_col = row['cluster_with_suffix']
		
		# Build list of subtypes for each cluster
		if cluster not in cluster_to_subtypes:
			cluster_to_subtypes[cluster] = []
			cluster_to_column[cluster] = cluster_col
			
		if subtype not in cluster_to_subtypes[cluster]:
			cluster_to_subtypes[cluster].append(subtype)
			
	return cluster_to_subtypes, cluster_to_column


def create_cohort_cluster_heatmaps(
	risk_contributions_per_subtype: pd.DataFrame,
	protective_contributions_per_subtype: pd.DataFrame,
	cluster_contributions_per_cohort: pd.DataFrame,
	cluster_subtype_mapping: pd.DataFrame,
	mode: str = 'combined',
	figsize_per_plot: Tuple[int, int] = (12, 5),
	risk_cmap: str = 'Reds',
	protective_cmap: str = 'Blues',
	combined_cmap: str = 'RdBu_r',
	save_path: str = None
) -> plt.Figure:
	"""
	Create heatmaps showing cohort contributions (risk and protective) to functional clusters 
	for each diabetes subtype.
	
	Parameters:
	-----------
	risk_contributions_per_subtype : pd.DataFrame
		DataFrame with subtypes as index and cohorts as columns.
		Values represent the risk weight of each cohort in each subtype.
	
	protective_contributions_per_subtype : pd.DataFrame
		DataFrame with subtypes as index and cohorts as columns.
		Values represent the protective weight of each cohort in each subtype.
	
	cluster_contributions_per_cohort : pd.DataFrame
		DataFrame with cohorts as index and functional clusters as columns.
		Column names should match the 'functional_cluster' values from cluster_subtype_mapping.
		Values represent the contribution of each cluster from each cohort.
	
	cluster_subtype_mapping : pd.DataFrame
		DataFrame with 'functional_cluster', 'subtype', and 'cluster_with_suffix' columns.
		Each row represents one cluster-subtype pairing.
		The 'functional_cluster' values should match column names in cluster_contributions_per_cohort.
		The 'cluster_with_suffix' is used for reference/tracking but not for data lookup.
		Example:
			functional_cluster | subtype | cluster_with_suffix
			beta_cell_expression | SIDD | beta_cell_expression_count
			obesity | MOD | obesity_count
			obesity | SIRD | obesity_count
	
	mode : str, optional
		'combined' - Show risk (positive) and protective (negative) in one heatmap per subtype
		'separate' - Show risk and protective as separate heatmaps side-by-side
		(default: 'combined')
	
	figsize_per_plot : Tuple[int, int], optional
		Size of each individual heatmap (default: (12, 5))
	
	risk_cmap : str, optional
		Colormap for risk contributions in separate mode (default: 'Reds')
	
	protective_cmap : str, optional
		Colormap for protective contributions in separate mode (default: 'Blues')
	
	combined_cmap : str, optional
		Colormap for combined mode (default: 'RdBu_r')
	
	save_path : str, optional
		Path to save the figure. If None, figure is not saved.
	
	Returns:
	--------
	fig : matplotlib.figure.Figure
		Figure object containing all heatmaps
	"""
	
	# Parse cluster-subtype mapping
	cluster_mapping, cluster_columns = parse_cluster_subtype_mapping(cluster_subtype_mapping)
	
	subtypes = risk_contributions_per_subtype.index.tolist()
	cohorts = risk_contributions_per_subtype.columns.tolist()
	n_subtypes = len(subtypes)
	
	if mode == 'combined':
		# Create subplots - one per subtype
		fig, axes = plt.subplots(n_subtypes, 1, 
								figsize=(figsize_per_plot[0], figsize_per_plot[1] * n_subtypes))
		if n_subtypes == 1:
			axes = [axes]
			
		fig.suptitle('Combined Risk (Positive) and Protective (Negative) Contributions\nto Functional Clusters by Diabetes Subtype', 
					fontsize=16, fontweight='bold', y=0.995)
		
		for idx, subtype in enumerate(subtypes):
			ax = axes[idx]
			_plot_combined_heatmap(
				subtype, cohorts, cluster_mapping, cluster_columns,
				risk_contributions_per_subtype, protective_contributions_per_subtype,
				cluster_contributions_per_cohort, ax, combined_cmap
			)
			
	else:  # mode == 'separate'
		# Create subplots - two columns (risk and protective) per subtype
		fig, axes = plt.subplots(n_subtypes, 2, 
								figsize=(figsize_per_plot[0] * 2, figsize_per_plot[1] * n_subtypes))
		if n_subtypes == 1:
			axes = axes.reshape(1, -1)
			
		fig.suptitle('Risk and Protective Contributions to Functional Clusters by Diabetes Subtype', 
					fontsize=16, fontweight='bold', y=0.995)
	
		for idx, subtype in enumerate(subtypes):
			ax_risk = axes[idx, 0]
			ax_protective = axes[idx, 1]
			
			_plot_risk_protective_heatmaps(
				subtype, cohorts, cluster_mapping, cluster_columns,
				risk_contributions_per_subtype, protective_contributions_per_subtype,
				cluster_contributions_per_cohort, ax_risk, ax_protective,
				risk_cmap, protective_cmap
			)
			
	plt.tight_layout()
	
	if save_path:
		plt.savefig(save_path, dpi=300, bbox_inches='tight')
		print(f"Figure saved to {save_path}")
		
	return fig


def _plot_combined_heatmap(subtype, cohorts, cluster_mapping, cluster_columns,
							risk_df, protective_df, cluster_df, ax, cmap):
	"""Helper function to plot combined risk/protective heatmap."""
	
	relevant_clusters = [
		cluster for cluster, subtypes_list in cluster_mapping.items()
		if subtype in subtypes_list
	]
	
	if not relevant_clusters:
		ax.text(0.5, 0.5, f'No functional clusters mapped to {subtype}', 
				ha='center', va='center', fontsize=12)
		ax.set_title(f'{subtype}', fontsize=14, fontweight='bold', pad=20)
		ax.axis('off')
		return
	
	# Calculate net contributions (risk - protective)
	contribution_matrix = np.zeros((len(cohorts), len(relevant_clusters)))
	
	for i, cohort in enumerate(cohorts):
		risk_weight = risk_df.loc[subtype, cohort]
		protective_weight = protective_df.loc[subtype, cohort]
		net_weight = risk_weight - protective_weight
		
		for j, cluster in enumerate(relevant_clusters):
			# Use the functional_cluster name directly (not the suffix version)
			cluster_contrib = cluster_df.loc[cohort, cluster]
			contribution_matrix[i, j] = net_weight * cluster_contrib
			
	# Create DataFrame for heatmap
	heatmap_df = pd.DataFrame(
		contribution_matrix,
		index=cohorts,
		columns=[c.replace('_', ' ').title() for c in relevant_clusters]
	)
	
	# Determine color scale limits (symmetric around 0)
	vmax = np.abs(contribution_matrix).max()
	vmin = -vmax if vmax > 0 else -0.001
	vmax = vmax if vmax > 0 else 0.001
	
	# Create heatmap
	sns.heatmap(
		heatmap_df,
		annot=True,
		fmt='.4f',
		cmap=cmap,
		center=0,
		vmin=vmin,
		vmax=vmax,
		cbar_kws={'label': 'Net Contribution (Risk - Protective)'},
		ax=ax,
		linewidths=0.5,
		linecolor='gray'
	)
	
	ax.set_title(f'{subtype}', fontsize=14, fontweight='bold', pad=20)
	ax.set_xlabel('Functional Cluster', fontsize=11, fontweight='bold')
	ax.set_ylabel('Cohort', fontsize=11, fontweight='bold')
	ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
	
	
def _plot_risk_protective_heatmaps(subtype, cohorts, cluster_mapping, cluster_columns,
									risk_df, protective_df, cluster_df,
									ax_risk, ax_protective, risk_cmap, protective_cmap):
	"""Helper function to plot separate risk and protective heatmaps."""
	
	relevant_clusters = [
		cluster for cluster, subtypes_list in cluster_mapping.items()
		if subtype in subtypes_list
	]
	
	if not relevant_clusters:
		for ax in [ax_risk, ax_protective]:
			ax.text(0.5, 0.5, f'No functional clusters mapped to {subtype}', 
					ha='center', va='center', fontsize=12)
			ax.axis('off')
		ax_risk.set_title(f'{subtype} - Risk', fontsize=14, fontweight='bold', pad=20)
		ax_protective.set_title(f'{subtype} - Protective', fontsize=14, fontweight='bold', pad=20)
		return
	
	# Calculate risk contributions
	risk_matrix = np.zeros((len(cohorts), len(relevant_clusters)))
	protective_matrix = np.zeros((len(cohorts), len(relevant_clusters)))
	
	for i, cohort in enumerate(cohorts):
		risk_weight = risk_df.loc[subtype, cohort]
		protective_weight = protective_df.loc[subtype, cohort]
		
		for j, cluster in enumerate(relevant_clusters):
			# Use the functional_cluster name directly (not the suffix version)
			cluster_contrib = cluster_df.loc[cohort, cluster]
			risk_matrix[i, j] = risk_weight * cluster_contrib
			protective_matrix[i, j] = protective_weight * cluster_contrib
			
	# Create DataFrames for heatmaps
	cluster_labels = [c.replace('_', ' ').title() for c in relevant_clusters]
	risk_df_heatmap = pd.DataFrame(risk_matrix, index=cohorts, columns=cluster_labels)
	protective_df_heatmap = pd.DataFrame(protective_matrix, index=cohorts, columns=cluster_labels)
	
	# Plot risk heatmap
	vmax_risk = risk_matrix.max() if risk_matrix.max() > 0 else 0.001
	sns.heatmap(
		risk_df_heatmap,
		annot=True,
		fmt='.4f',
		cmap=risk_cmap,
		cbar_kws={'label': 'Risk Contribution'},
		ax=ax_risk,
		linewidths=0.5,
		linecolor='gray',
		vmin=0,
		vmax=vmax_risk
	)
	
	ax_risk.set_title(f'{subtype} - Risk', fontsize=14, fontweight='bold', pad=20)
	ax_risk.set_xlabel('Functional Cluster', fontsize=11, fontweight='bold')
	ax_risk.set_ylabel('Cohort', fontsize=11, fontweight='bold')
	ax_risk.set_xticklabels(ax_risk.get_xticklabels(), rotation=45, ha='right')
	ax_risk.set_yticklabels(ax_risk.get_yticklabels(), rotation=0)
	
	# Plot protective heatmap
	vmax_protective = protective_matrix.max() if protective_matrix.max() > 0 else 0.001
	sns.heatmap(
		protective_df_heatmap,
		annot=True,
		fmt='.4f',
		cmap=protective_cmap,
		cbar_kws={'label': 'Protective Contribution'},
		ax=ax_protective,
		linewidths=0.5,
		linecolor='gray',
		vmin=0,
		vmax=vmax_protective
	)
	
	ax_protective.set_title(f'{subtype} - Protective', fontsize=14, fontweight='bold', pad=20)
	ax_protective.set_xlabel('Functional Cluster', fontsize=11, fontweight='bold')
	ax_protective.set_ylabel('Cohort', fontsize=11, fontweight='bold')
	ax_protective.set_xticklabels(ax_protective.get_xticklabels(), rotation=45, ha='right')
	ax_protective.set_yticklabels(ax_protective.get_yticklabels(), rotation=0)
	
	

if __name__ == "__main__":
	"""
	Complete example of how to use all subtype enrichment functions
	"""
	
	print("=" * 80)
	print("SUBTYPE ENRICHMENT ANALYSIS - COMPLETE WORKFLOW")
	print("=" * 80)
	
	pheno_path = '/Users/kerimulterer/prsInteractive/results/type2Diabetes'
	# Load your data
	df = pd.read_csv(f'{pheno_path}/scores/clusterWeights.csv')  # or pd.read_csv() if it's actually CSV
	
	# Run the analysis
	results = rank_functional_clusters_by_cohort(df)
	
	# Display results
#	display_results(results)
	display_direction_analysis(results)
	
	signatures = pd.DataFrame()
	for cohort in results['scores_df'].index:
		signature = find_cohort_signature_clusters(results,cohort,top_n=results['scores_df'].shape[1])
		signature['cohort'] = cohort
		signatures = pd.concat([signature,signatures],ignore_index=True)
	

			


	# Step 1: Load cluster-subtype mapping
	print("\nLoading cluster-subtype mapping...")
	cluster_mapping = load_cluster_subtype_mapping(f'{pheno_path}/scores/clusterToSubtype.csv')
	
	print("=" * 80)
	print("COMPARING WEIGHTING METHODS")
	print("=" * 80)
	
	# Method 1: Effect size weighting
	print("\n1. Effect Size Weighting (recommended)")
	subtype_weighted = calculate_weighted_subtype_enrichment(
		df, results, cluster_mapping,
		weighting='effect_size'
	)
	display_weighted_results(subtype_weighted)
	plot_weighted_subtype_pies(subtype_weighted, f'{pheno_path}/figures/weighted_by_effect_size.png')
	
	# Method 2: Distinctiveness weighting
	print("\n2. Distinctiveness Weighting (effect × uniqueness)")
	subtype_distinctive = calculate_weighted_subtype_enrichment(
		df, results, cluster_mapping,
		weighting='distinctiveness'
	)
	plot_weighted_subtype_pies(subtype_distinctive, f'{pheno_path}/figures/weighted_by_distinctiveness.png')
	
	# Method 3: Normalized weighting
	print("\n3. Normalized Score Weighting")
	subtype_normalized = calculate_weighted_subtype_enrichment(
		df, results, cluster_mapping,
		weighting='normalized'
	)
	plot_weighted_subtype_pies(subtype_normalized, f'{pheno_path}/figures/weighted_by_normalized.png')
	
	# Create comparison
	print("\n4. Creating Comparison Plots...")
	fig, comparison = plot_contribution_comparison(subtype_weighted)


	# Export results to Excel with multiple sheets
	with pd.ExcelWriter(f'{pheno_path}/scores/cluster_rankings_results.xlsx') as writer:
		results['scores_df'].to_excel(writer, sheet_name='Normalized Scores')
		results['cluster_ranks_df'].to_excel(writer, sheet_name='Cluster Rankings')
		results['cohort_ranks_df'].to_excel(writer, sheet_name='Cohort Rankings')
		results['top_cohorts'].to_excel(writer, sheet_name='Top Cohorts by Cluster')
		signatures.to_excel(writer,sheet_name='Signature Clusters by Cohort')
		
	print("\nResults exported to 'cluster_rankings_results.xlsx'")
	
		
	with pd.ExcelWriter(f'{pheno_path}/scores/cluster_subtype_results_distinctiveness.xlsx') as writer:
		subtype_distinctive['cohort_weights_per_subtype'].to_excel(writer, sheet_name='cohort_weights_per_subtype')
		subtype_distinctive['cohort_percentages_per_subtype'].to_excel(writer, sheet_name='cohort_percentages_per_subtype')
		subtype_distinctive['gene_counts'].to_excel(writer, sheet_name='gene_counts')
		subtype_distinctive['risk_contributions'].to_excel(writer, sheet_name='risk_contributions')
		subtype_distinctive['protective_contributions'].to_excel(writer, sheet_name='protective_contributions')
		
	print("\nResults exported to 'cluster_subtype_results_distinctiveness.xlsx'")
	
		
	with pd.ExcelWriter(f'{pheno_path}/scores/cluster_subtype_results_effect_sizes.xlsx') as writer:
		subtype_weighted['cohort_weights_per_subtype'].to_excel(writer, sheet_name='cohort_weights_per_subtype')
		subtype_weighted['cohort_percentages_per_subtype'].to_excel(writer, sheet_name='cohort_percentages_per_subtype')
		subtype_weighted['gene_counts'].to_excel(writer, sheet_name='gene_counts')
		subtype_weighted['risk_contributions'].to_excel(writer, sheet_name='risk_contributions')
		subtype_weighted['protective_contributions'].to_excel(writer, sheet_name='protective_contributions')
		
	print("\nResults exported to 'cluster_subtype_results_effect_sizes.xlsx'")
	
	with pd.ExcelWriter(f'{pheno_path}/scores/cluster_subtype_results_normalized_gene_count.xlsx') as writer:
		subtype_normalized['cohort_weights_per_subtype'].to_excel(writer, sheet_name='cohort_weights_per_subtype')
		subtype_normalized['cohort_percentages_per_subtype'].to_excel(writer, sheet_name='cohort_percentages_per_subtype')
		subtype_normalized['gene_counts'].to_excel(writer, sheet_name='gene_counts')
		subtype_normalized['risk_contributions'].to_excel(writer, sheet_name='risk_contributions')
		subtype_normalized['protective_contributions'].to_excel(writer, sheet_name='protective_contributions')
	
	print("\nResults exported to 'cluster_subtype_results_normalized_gene_count.xlsx'")
	


	
	print("\n5. Creating composition and risk_protective comparison Plots...")
	
	scores_df_plot = results['scores_df']
	scores_df_plot.columns = [col.replace("_count","") for col in scores_df_plot.columns]
	
	#plot heatmaps
	fig = create_cohort_cluster_heatmaps(
		subtype_weighted['risk_contributions'],
		subtype_weighted['protective_contributions'],
		results['scores_df'],
		cluster_mapping,
		'combined',
		(12, 5),
		'Reds',
		'Blues',
		'RdBu_r',
		save_path=f'{pheno_path}/figures/cohort_composition_heatmap_combined.png'
	)
	
	#plot heatmaps
	fig = create_cohort_cluster_heatmaps(
		subtype_weighted['risk_contributions'],
		subtype_weighted['protective_contributions'],
		results['scores_df'],
		cluster_mapping,
		'separate',
		(12, 5),
		'Reds',
		'Blues',
		'RdBu_r',
		save_path=f'{pheno_path}/figures/cohort_composition_heatmap_separated.png'
	)
		
#	plot_signature_distinctiveness_scatter(results,save_path=f'{pheno_path}/figures/subtype_composition_by_cohort.png')
	plot_risk_protective_comparison(results, save_path=f'{pheno_path}/figures/risk_protective_comparison.png')
	
	print("\n" + "=" * 80)
	print("ANALYSIS COMPLETE!")
	print("=" * 80)
	print("\nOutputs created:")
	print("  ✓ Console: Detailed subtype analysis")
	print("  ✓ subtype_composition_by_cohort.png: 5 pie charts (custom colors)")
	print("  ✓ weighted_by_effect_size.png")
	print("  ✓ weighted_by_distinctiveness.png")
	print("  ✓ weighted_by_normalized.png")
	print("  ✓ contribution_comparison.png")
	print("  ✓ signature_distinctiveness_scatter.png: Scatter plot of top clusters")
	print("  ✓ risk_protective_comparison.png: bar chart of risk protection of top cluster comparisons")