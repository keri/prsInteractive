#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from typing import Dict, List, Tuple, Optional
from matplotlib.colors import LinearSegmentedColormap


# ============================================================================
# CALCULATE FUNCTIONAL CLUSTER RANKING ENRICHMENT
# ============================================================================

# Colorblind-friendly palettes
# Okabe-Ito palette for cohorts
COHORT_COLORS = {
	'G': '#E69F00',      # Orange
	'GxG': '#56B4E9',    # Sky blue
	'GxGxE': '#009E73'   # Bluish green
}

# Colorblind-friendly palette for subtypes (using ColorBrewer Set2)
SUBTYPE_COLORS = {
	'SIRD': '#8DA0CB',   # Light blue-purple
	'MOD': '#FC8D62',    # Salmon
	'SIDD': '#66C2A5',   # Teal
	'MARD': '#E78AC3',   # Pink
	'SAID': '#A6D854'    # Yellow-green
}

# Cohort markers for additional distinction
COHORT_MARKERS = {
	'G': 'o',      # Circle
	'GxG': 's',    # Square
	'GxGxE': '^'   # Triangle
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

	# Calculate percentages based on RISK contributions only
	percentages_df = risk_df.div(risk_df.sum(axis=1), axis=0) * 100
	percentages_df = percentages_df.fillna(0)  # Handle cases where risk sum is 0

	# Summary
	enrichment_summary = {
			'total_subtypes': len(risk_df),
			'cohorts_contributing': risk_df.columns.tolist(),
			'weighting_method': weighting,
			'total_weighted_contribution_per_subtype': risk_df.sum(axis=1).to_dict()
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

	# Colorblind-friendly palettes
	# Okabe-Ito palette for cohorts
	COHORT_COLORS = {
		'G': '#E69F00',      # Orange
		'GxG': '#56B4E9',    # Sky blue
		'GxGxE': '#009E73'   # Bluish green
	}
	
	# Colorblind-friendly palette for subtypes (using ColorBrewer Set2)
	SUBTYPE_COLORS = {
		'SIRD': '#8DA0CB',   # Light blue-purple
		'MOD': '#FC8D62',    # Salmon
		'SIDD': '#66C2A5',   # Teal
		'MARD': '#E78AC3',   # Pink
		'SAID': '#A6D854'    # Yellow-green
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

#def calculate_subtype_enrichment(df, results, cluster_subtype_mapping, 
#									cohort_col='cohort'):
#	"""
#	Calculate T2D subtype enrichment showing cohort composition within each subtype.
#	
#	This determines what proportion of each subtype's associations come from each cohort.
#	
#	Parameters:
#	-----------
#	df : pandas.DataFrame
#		Original data with genes and cluster associations
#	results : dict
#		Results from rank_functional_clusters_by_cohort
#	cluster_subtype_mapping : DataFrame
#		Mapping from load_cluster_subtype_mapping()
#	cohort_col : str
#		Name of cohort column
#		
#	Returns:
#	--------
#	dict with:
#		- 'cohort_counts_per_subtype': DataFrame (subtypes × cohorts) with counts
#		- 'cohort_percentages_per_subtype': DataFrame (subtypes × cohorts) with percentages
#		- 'enrichment_summary': Summary statistics
#	
#	Example:
#	--------
#	>>> mapping = load_cluster_subtype_mapping('clusterToSubtype.csv')
#	>>> subtype_results = calculate_subtype_enrichment(df, results, mapping)
#	>>> print(subtype_results['cohort_percentages_per_subtype'])
#	"""
#	
#	
#	df_expanded,cluster_cols = expand_cohort_rows(df, 'cohort')
#	
#	# Create mapping dictionary for quick lookup
#	cluster_to_subtypes = {}
#	for cluster in cluster_cols:
#		cluster_name = cluster.replace('_count', '')
#		matching = cluster_subtype_mapping[
#			cluster_subtype_mapping['functional_cluster'] == cluster_name
#		]
#		if len(matching) > 0:
#			cluster_to_subtypes[cluster] = matching['subtype'].tolist()
#		else:
#			cluster_to_subtypes[cluster] = []
#			
#	# Calculate subtype-cohort counts (what cohorts contribute to each subtype)
#	subtype_cohort_counts = {}
#	
#	for cohort in df_expanded[cohort_col].unique():
#		cohort_data = df_expanded[df_expanded[cohort_col] == cohort]
#		
#		# For each gene in the cohort
#		for _, row in cohort_data.iterrows():
#			# Check each cluster
#			for cluster in cluster_cols:
#				# If gene is associated with this cluster
#				if pd.notna(row.get(cluster)) and row[cluster] == 1:
#					# Get subtypes for this cluster
#					subtypes = cluster_to_subtypes.get(cluster, [])
#					
#					# Count this cohort's contribution to each subtype
#					for subtype in subtypes:
#						if subtype not in subtype_cohort_counts:
#							subtype_cohort_counts[subtype] = {}
#						if cohort not in subtype_cohort_counts[subtype]:
#							subtype_cohort_counts[subtype][cohort] = 0
#						subtype_cohort_counts[subtype][cohort] += 1
#						
#	# Convert to DataFrame (rows = subtypes, columns = cohorts)
#	counts_df = pd.DataFrame(subtype_cohort_counts).T.fillna(0).astype(int)
#	counts_df.index.name = 'subtype'
#	
#	# Calculate percentages (what % of each subtype comes from each cohort)
#	percentages_df = counts_df.div(counts_df.sum(axis=1), axis=0) * 100
#	
#	# Summary statistics
#	enrichment_summary = {
#		'total_subtypes': len(counts_df),
#		'cohorts_contributing': counts_df.columns.tolist(),
#		'total_associations_per_subtype': counts_df.sum(axis=1).to_dict()
#	}
#	
#	return {
#		'cohort_counts_per_subtype': counts_df,
#		'cohort_percentages_per_subtype': percentages_df,
#		'enrichment_summary': enrichment_summary
#	}
	

	

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

#def display_subtype_results(subtype_results):
#	"""
#	Display subtype enrichment results in readable format.
#	
#	Parameters:
#	-----------
#	subtype_results : dict
#		Results from calculate_subtype_enrichment
#	
#	Example:
#	--------
#	>>> display_subtype_results(subtype_results)
#	"""
#	
#	print("\n" + "=" * 80)
#	print("T2D SUBTYPE ENRICHMENT ANALYSIS")
#	print("=" * 80)
#	
#	print(f"\nSubtypes identified: {', '.join(subtype_results['cohort_counts_per_subtype'].index.tolist())}")
#	print(f"Total cohorts analyzed: {len(subtype_results['cohort_counts_per_subtype'].columns)}")
#	
#	print("\n" + "=" * 80)
#	print("COHORT COUNTS WITHIN EACH SUBTYPE")
#	print("(How many associations from each cohort contribute to each subtype)")
#	print("=" * 80)
#	print(subtype_results['cohort_counts_per_subtype'].to_string())
#	
#	print("\n" + "=" * 80)
#	print("COHORT PERCENTAGES WITHIN EACH SUBTYPE")
#	print("(What % of each subtype comes from each cohort)")
#	print("=" * 80)
#	print(subtype_results['cohort_percentages_per_subtype'].round(1).to_string())
#	
#	# Find dominant cohort per subtype
#	print("\n" + "=" * 80)
#	print("DOMINANT COHORT FOR EACH SUBTYPE")
#	print("=" * 80)
#	
#	for subtype, row in subtype_results['cohort_percentages_per_subtype'].iterrows():
#		if row.sum() > 0:
#			dominant = row.idxmax()
#			percentage = row.max()
#			total = subtype_results['cohort_counts_per_subtype'].loc[subtype].sum()
#			print(f"{subtype}: {dominant} contributes {percentage:.1f}% (n={int(total)} total associations)")
#			
#	# Show subtype composition details
#	print("\n" + "=" * 80)
#	print("SUBTYPE COMPOSITION DETAILS")
#	print("=" * 80)
#	
#	for subtype in subtype_results['cohort_percentages_per_subtype'].index:
#		total = subtype_results['cohort_counts_per_subtype'].loc[subtype].sum()
#		print(f"\n{subtype} (n={int(total)} associations):")
#		
#		composition = subtype_results['cohort_percentages_per_subtype'].loc[subtype].sort_values(ascending=False)
#		for cohort, pct in composition.items():
#			if pct > 0:
#				count = subtype_results['cohort_counts_per_subtype'].loc[subtype, cohort]
#				print(f"  {cohort}: {pct:.1f}% (n={int(count)})")
#
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
	
	
def prepare_forest_plot_data(risk_contributions_per_subtype,
	protective_contributions_per_subtype,
	cluster_contributions_per_cohort,
	cluster_subtype_mapping):
	"""
		Prepare data for forest plot showing functional cluster contributions by cohort for each subtype.
		
		Parameters:
		-----------
		risk_contributions_per_subtype : pd.DataFrame
				Risk weights per cohort per subtype
		protective_contributions_per_subtype : pd.DataFrame
				Protective weights per cohort per subtype
		cluster_contributions_per_cohort : pd.DataFrame
				Cohorts as index, functional clusters as columns
		cluster_subtype_mapping : pd.DataFrame
				Columns: 'functional_cluster', 'subtype', 'cluster_with_suffix'
		
		Returns:
		--------
		dict : {subtype: {cluster: {cohort: net_contribution}}}
		"""
	forest_data = {}
	
	# Get net contributions (risk - protective) for weighting
	net_cohort_contributions = risk_contributions_per_subtype 
#	net_cohort_contributions = risk_contributions_per_subtype - protective_contributions_per_subtype
	
	cluster_contributions_per_cohortCopy = cluster_contributions_per_cohort.copy()
	cluster_contributions_per_cohortCopy.columns = [col.replace('_count','') for col in cluster_contributions_per_cohort.columns]
	
	# Group mapping by subtype
	for subtype in cluster_subtype_mapping['subtype'].unique():
		forest_data[subtype] = {}
		
		# Get clusters for this subtype
		subtype_clusters = cluster_subtype_mapping[
			cluster_subtype_mapping['subtype'] == subtype
		]['functional_cluster'].unique()
		
		for cluster in subtype_clusters:
			# Skip if cluster not in contributions dataframe
			if cluster not in cluster_contributions_per_cohortCopy.columns:
				continue
			
			forest_data[subtype][cluster] = {}
			
			# For each cohort, calculate weighted contribution
			for cohort in cluster_contributions_per_cohortCopy.index:
				# Get cluster contribution from this cohort
				cluster_contrib = cluster_contributions_per_cohortCopy.loc[cohort, cluster]
				
				# Weight by the cohort's net contribution to this subtype
				if cohort in net_cohort_contributions.columns and subtype in net_cohort_contributions.index:
					cohort_weight = net_cohort_contributions.loc[subtype, cohort]
					weighted_contrib = cluster_contrib * cohort_weight
				else:
					weighted_contrib = 0
					
				forest_data[subtype][cluster][cohort] = weighted_contrib
						
	return forest_data
	
#def prepare_forest_plot_data(risk_contributions_per_subtype,
#	protective_contributions_per_subtype,
#	cluster_contributions_per_cohort,
#	cluster_subtype_mapping):
#	"""
#	Prepare data for forest plot showing functional cluster contributions by cohort for each subtype.
#	Now preserves separate risk and protective contributions.
#	
#	Parameters:
#	-----------
#	risk_contributions_per_subtype : pd.DataFrame
#			Risk weights per cohort per subtype
#	protective_contributions_per_subtype : pd.DataFrame
#			Protective weights per cohort per subtype
#	cluster_contributions_per_cohort : pd.DataFrame
#			Cohorts as index, functional clusters as columns
#	cluster_subtype_mapping : pd.DataFrame
#			Columns: 'functional_cluster', 'subtype', 'cluster_with_suffix'
#	
#	Returns:
#	--------
#	dict : {subtype: {cluster: {cohort: {'risk': value, 'protective': value}}}}
#	"""
#	forest_data = {}
#	
#	cluster_contributions_per_cohortCopy = cluster_contributions_per_cohort.copy()
#	cluster_contributions_per_cohortCopy.columns = [col.replace('_count','') for col in cluster_contributions_per_cohort.columns]
#	
#	# Group mapping by subtype
#	for subtype in cluster_subtype_mapping['subtype'].unique():
#		forest_data[subtype] = {}
#		
#		# Get clusters for this subtype
#		subtype_clusters = cluster_subtype_mapping[
#			cluster_subtype_mapping['subtype'] == subtype
#		]['functional_cluster'].unique()
#		
#		for cluster in subtype_clusters:
#			# Skip if cluster not in contributions dataframe
#			if cluster not in cluster_contributions_per_cohortCopy.columns:
#				continue
#			
#			forest_data[subtype][cluster] = {}
#			
#			# For each cohort, calculate weighted contributions (keep risk and protective separate)
#			for cohort in cluster_contributions_per_cohortCopy.index:
#				# Get cluster contribution from this cohort
#				cluster_contrib = cluster_contributions_per_cohortCopy.loc[cohort, cluster]
#				
#				# Weight by the cohort's risk and protective contributions to this subtype
#				risk_weighted = 0
#				protective_weighted = 0
#				
#				if cohort in risk_contributions_per_subtype.columns and subtype in risk_contributions_per_subtype.index:
#					risk_weight = risk_contributions_per_subtype.loc[subtype, cohort]
#					risk_weighted = cluster_contrib * risk_weight
#					
#					if cohort in protective_contributions_per_subtype.columns and subtype in protective_contributions_per_subtype.index:
#						protective_weight = protective_contributions_per_subtype.loc[subtype, cohort]
#						protective_weighted = cluster_contrib * protective_weight
#						
#						forest_data[subtype][cluster][cohort] = {
#							'risk': risk_weighted,
#							'protective': protective_weighted
#						}
#						
#	return forest_data
#	
	
#def create_forest_plot(forest_data, output_path, top_n_clusters=10):
#	"""Create forest plot showing functional cluster contributions"""
#	
#	fig, ax_forest = plt.subplots(figsize=(14, 16))
#	
#	subtypes = ['SIRD', 'MOD', 'SIDD', 'MARD', 'SAID']
#	cohorts = ['G', 'GxG', 'GxGxE']
#	
#	# Prepare data for forest plot
#	y_positions = []
#	x_values = []
#	colors_list = []
#	markers_list = []
#	
#	y_pos = 0
#	y_ticks = []
#	y_labels = []
#	subtype_positions = {}  # Track subtype y-positions for horizontal lines
#	
#	
#	# Iterate through subtypes
#	for subtype in subtypes:
#		
#		# Store subtype position for horizontal line
#		subtype_positions[subtype] = y_pos
#		y_pos -= 0.8
#		
#		# Add subtype header
#		y_ticks.append(y_pos)
#		y_labels.append(f'{subtype}')
#		y_pos -= 0.8
#		
#		# Get clusters for this subtype
#		if subtype not in forest_data or not forest_data[subtype]:
#			y_pos -= 0.5
#			continue
#		
#		clusters = forest_data[subtype]
#		
#		# Sort clusters by total absolute contribution (for better visualization)
#		cluster_totals = {
#			cluster: sum(abs(forest_data.get(cohort, 0)['risk'] - forest_data.get(cohort, 0)['protective']) 
#				for cohort in cohorts)
#			for cluster, forest_data in clusters.items()
#		}
#		sorted_clusters = sorted(cluster_totals.items(), 
#			key=lambda x: x[1], 
#			reverse=True)
#		
#		# Plot each cluster
#		for cluster, total_contrib in sorted_clusters[:top_n_clusters]:
#			if total_contrib < 1e-6:  # Skip clusters with negligible contribution
#				continue
#		
#		cohort_data = clusters[cluster]
#		has_contribution = False
#		
#		# Plot each cohort's contribution for this cluster
#		for cohort in cohorts:
#			contribution = abs(cohort_data.get(cohort, 0)['risk'] - cohort_data.get(cohort, 0)['protective'])
#			
#			if abs(contribution) > 0.0001:  # Only plot non-zero contributions
#				y_positions.append(y_pos)
#				x_values.append(contribution)
#				colors_list.append(COHORT_COLORS[cohort])
#				markers_list.append(COHORT_MARKERS[cohort])
#				has_contribution = True
#			
#			# Only add cluster label if it has contributions
#			if has_contribution:
#				y_ticks.append(y_pos)
#				# Clean up cluster name
#				cluster_display = cluster.replace('_', ' ').title()
#				if len(cluster_display) > 35:
#					cluster_display = cluster_display[:32] + '...'
#					y_labels.append(f'  {cluster_display}')
#					y_pos -= 1.0
#					
#					y_pos -= 0.5  # Extra space between subtypes
#					
#		# Get x-axis limits first (before plotting) to span full width
#		if x_values:
#			x_min = min(x_values)
#			x_max = max(x_values)
#			x_range = x_max - x_min
#			x_padding = x_range * 0.1
#			xlim = (x_min - x_padding, x_max + x_padding)
#			ax_forest.set_xlim(xlim)
#		else:
#			xlim = ax_forest.get_xlim()
#			
#			# Draw horizontal lines for each subtype section
#		for subtype, y_subtype in subtype_positions.items():
#			# Draw horizontal line spanning the plot
#			ax_forest.axhline(y=y_subtype, color='gray', linestyle='-', 
#				linewidth=2, alpha=0.4, zorder=0)
#			
#			# Add subtype label centered on the line
#			x_center = (xlim[0] + xlim[1]) / 2
#			ax_forest.text(x_center, y_subtype, f'  {subtype}  ',
#				ha='center', va='center',
#				fontsize=16, weight='bold',
#				color=SUBTYPE_COLORS[subtype],
#				bbox=dict(boxstyle='round,pad=0.5', 
#					facecolor='white', 
#					edgecolor=SUBTYPE_COLORS[subtype],
#					linewidth=2,
#					alpha=0.95),
#				zorder=5)
#				
#				# Add horizontal lines connecting to dots
#			for y, x in zip(y_positions, x_values):
#				ax_forest.plot([0, x], [y, y], 
#					color='gray', 
#					linestyle='-', 
#					linewidth=1.2, 
#					alpha=0.4,
#					zorder=1)
#				
#				# Plot the forest plot points
#			for i, (y, x, color, marker) in enumerate(zip(y_positions, x_values, colors_list, markers_list)):
#				ax_forest.scatter(x, y, 
#					marker=marker,
#					s=150, 
#					color=color,
#					edgecolor='black', 
#					linewidth=1.0,
#					alpha=0.85,
#					zorder=3)
#				
#				# Add vertical line at x=0
#			ax_forest.axvline(x=0, color='black', linestyle='-', linewidth=2, 
#				zorder=1, alpha=0.9)
#			
#			# Styling
#			ax_forest.set_yticks(y_ticks)
#			ax_forest.set_yticklabels(y_labels, fontsize=11)
#			ax_forest.set_xlabel('Net Contribution (Weighted)', 
#				fontsize=14, weight='bold')
#			ax_forest.set_title('Functional Cluster Contributions to T2D Subtypes',
#				fontsize=16, weight='bold', pad=20)
#			
#			# Style y-axis labels - bold for subtypes
#			for tick_label, label_text in zip(ax_forest.get_yticklabels(), y_labels):
#				if not label_text.startswith('  '):  # Subtype labels
#					tick_label.set_weight('bold')
#					tick_label.set_fontsize(16)
#					# Color code subtype labels
#					subtype_name = label_text.strip()
#				if subtype_name in SUBTYPE_COLORS:
#					tick_label.set_color(SUBTYPE_COLORS[subtype_name])
#					
#					# Grid
#			ax_forest.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.8, zorder=0)
#			ax_forest.set_axisbelow(True)
#			
#			# Spines
#			ax_forest.spines['top'].set_visible(False)
#			ax_forest.spines['right'].set_visible(False)
#			ax_forest.spines['left'].set_linewidth(1.5)
#			ax_forest.spines['bottom'].set_linewidth(1.5)
#			
#			# Add shaded regions to distinguish risk vs protective
#			#		xlim = ax_forest.get_xlim()
#			#		ax_forest.axvspan(0, xlim[1], alpha=0.05, color='white', zorder=0)
#			#		ax_forest.axvspan(xlim[0], 0, alpha=0.05, color='white', zorder=0)
#			
#			# Add text labels for risk/protective
#			y_max = y_ticks[0] if y_ticks else 0
#			ax_forest.text(xlim[1]*0.95, y_max + 1, 'Risk →', 
#				ha='right', va='bottom', fontsize=20, 
#				weight='bold', style='italic', color='darkred')
#			ax_forest.text(xlim[0]*0.95, y_max + 1, '← Protective', 
#				ha='left', va='bottom', fontsize=20, 
#				weight='bold', style='italic', color='darkblue')
#			
#			# Legend for cohorts
#			legend_elements = [
#				plt.Line2D([0], [0], marker=COHORT_MARKERS[cohort], 
#					color='w', 
#					markerfacecolor=COHORT_COLORS[cohort],
#					markeredgecolor='black',
#					markeredgewidth=1.0,
#					markersize=11, 
#					label=f'{cohort}',
#					linestyle='None')
#				for cohort in cohorts
#			]
#			
#			legend = ax_forest.legend(handles=legend_elements, 
#				loc='lower right',
#				title='ePRS model Type',
#				frameon=True,
#				fancybox=True,
#				shadow=True,
#				fontsize=11,
#				title_fontsize=12)
#			legend.get_frame().set_alpha(0.95)
#			legend.get_frame().set_linewidth(1.5)
#			
#			# Adjust layout
#			plt.tight_layout()
#			
#			# Save figure
#			plt.savefig(output_path, dpi=300, bbox_inches='tight', 
#				facecolor='white', edgecolor='none')
#			print(f"Forest plot saved: {output_path}")
#			plt.close()
				
	
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
	
def create_combined_visualization(df, forest_data, output_path, top_n_clusters=10):
		"""
		Create combined visualization with pie charts and forest plots for each subtype.
		
		Parameters:
		-----------
		df : dataFrame
				index : subtype,
				columns: cohorts,
				values: percentage contribution
		forest_data : dict
			{subtype: {cluster: {cohort: contribution}}}
		output_path : str
				Path to save the figure
		top_n_clusters : int
				Number of top clusters to show per subtype
		"""
	
		subtypes = df.index.tolist()
		cohorts = df.columns.tolist()
	
		# Create figure with GridSpec for better control
		fig = plt.figure(figsize=(18, 22))
		gs = gridspec.GridSpec(5, 2, figure=fig, width_ratios=[1, 2.5], 
													hspace=0.35, wspace=0.3,
													left=0.08, right=0.96, top=0.96, bottom=0.04)
	
		# Add main title
		fig.suptitle('Cohort Composition and Functional Cluster Contributions to T2D Subtypes',
								fontsize=20, weight='bold', y=0.985)
						
		# Convert to dict with proportions (divide by 100 since values are percentages)
		cohort_proportions = {}
		for subtype in df.index:
			cohort_proportions[subtype] = {
				cohort: df.loc[subtype, cohort]
				for cohort in df.columns
			}
	
		# Process each subtype
		for idx, subtype in enumerate(subtypes):
				# Create axes for this row
				ax_pie = fig.add_subplot(gs[idx, 0])
				ax_forest = fig.add_subplot(gs[idx, 1])
			
				# ===== PIE CHART =====
				values = [cohort_proportions[subtype].get(cohort, 0) for cohort in cohorts]
				colors = [COHORT_COLORS[cohort] for cohort in cohorts]
			
				# Only show labels for non-zero values
				labels = []
				for cohort, val in zip(cohorts, values):
						if val > 0.01:  # Show label if > 1%
								labels.append(cohort)
						else:
								labels.append('')
							
				wedges, texts = ax_pie.pie(
						values,
#						labels=labels,
						colors=colors
#						autopct=lambda pct: f'{pct:.1f}%' if pct > 1 else '',
#						startangle=90,
#						textprops={'fontsize': 11, 'weight': 'bold'}
				)
			
				# Make percentage text white for better visibility
#				for autotext in autotexts:
#						autotext.set_color('white')
#						autotext.set_fontsize(10)
					
#				ax_pie.set_title('Cohort Composition', fontsize=13, weight='bold', pad=10)
			
				# ===== FOREST PLOT =====
				if subtype not in forest_data or not forest_data[subtype]:
						ax_forest.text(0.5, 0.5, 'No data available', 
													ha='center', va='center', fontsize=12, style='italic')
						ax_forest.set_xlim(-1, 1)
						ax_forest.set_ylim(-1, 1)
						ax_forest.axis('off')
				else:
						clusters = forest_data[subtype]
					
					#Sort clusters by total absolute contribution
						cluster_totals = {
							cluster: sum(abs(cohort_data.get(cohort, 0)) 
								for cohort in cohorts)
							for cluster, cohort_data in clusters.items()
						}
					
#						cluster_totals = {
#							cluster: sum(abs(cohort_data.get(cohort, {}).get('risk', 0)) + 
#								abs(cohort_data.get(cohort, {}).get('protective', 0))
#								for cohort in cohorts)
#							for cluster, cohort_data in clusters.items()
#						}
						sorted_clusters = sorted(cluster_totals.items(), 
													key=lambda x: x[1], 
													reverse=True)
					
						# Prepare data for plotting
						y_positions = []
						x_values = []
						colors_list = []
						markers_list = []
						y_ticks = []
						y_labels = []
					
						y_pos = 0
					
						# Plot each cluster
						for cluster, total_contrib in sorted_clusters[:top_n_clusters]:
								if total_contrib < 1e-6:
										continue
							
								cohort_data = clusters[cluster]
								has_contribution = False
							
								# Plot each cohort's contribution
								for cohort in cohorts:
									contribution = cohort_data.get(cohort, 0)
										
									if cohort not in cohort_data:
										continue
									
#									risk_val = cohort_data[cohort].get('risk', 0)
#									protective_val = cohort_data[cohort].get('protective', 0)
									
									if abs(contribution) > 0.0001:
										y_positions.append(y_pos)
										x_values.append(contribution)
										colors_list.append(COHORT_COLORS[cohort])
										markers_list.append(COHORT_MARKERS[cohort])
										has_contribution = True
									
									# Plot risk contribution (positive, on right side)
#									if abs(risk_val) > 0.0001:
#										y_positions.append(y_pos)
#										x_values.append(risk_val)
#										colors_list.append(COHORT_COLORS[cohort])
#										markers_list.append(COHORT_MARKERS[cohort])
#										has_contribution = True
#										
#										# Plot protective contribution (negative, on left side)
#										if abs(protective_val) > 0.0001:
#											y_positions.append(y_pos)
#											x_values.append(-protective_val)  # Negative for left side
#											colors_list.append(COHORT_COLORS[cohort])
#											markers_list.append(COHORT_MARKERS[cohort])
#											has_contribution = True
											
								if has_contribution:
										y_ticks.append(y_pos)
										# Clean up cluster name
										cluster_display = cluster.replace('_', ' ')
										if len(cluster_display) > 35:
												cluster_display = cluster_display[:32] + '...'
										y_labels.append(cluster_display)
										y_pos -= 1.0
									
						# Get x-axis limits
						if x_values:
#								x_min = min(x_values)
								x_min = 0
								x_max = max(x_values)
								x_range = x_max - x_min
								if x_range > 0:
										x_padding = x_range * 0.15
#										xlim = (x_min - x_padding, x_max + x_padding)
										xlim = (x_min, x_max + x_padding)
								else:
										xlim = (-1, 1)
								ax_forest.set_xlim(xlim)
						else:
								xlim = (-1, 1)
								ax_forest.set_xlim(xlim)
							
						# Add horizontal lines connecting to dots
						for y, x in zip(y_positions, x_values):
								ax_forest.plot([0, x], [y, y], 
															color='gray', 
															linestyle='-', 
															linewidth=1.0, 
															alpha=0.3,
															zorder=1)
							
						# Plot the forest plot points
						for y, x, color, marker in zip(y_positions, x_values, colors_list, markers_list):
								ax_forest.scatter(x, y, 
													marker=marker,
													s=120, 
													color=color,
													edgecolor='black', 
													linewidth=0.8,
													alpha=0.85,
													zorder=3)
							
						# Add vertical line at x=0
						ax_forest.axvline(x=0, color='black', linestyle='-', linewidth=1.5, 
															zorder=2, alpha=0.8)
					
						# Styling
						if y_ticks:
								ax_forest.set_yticks(y_ticks)
								ax_forest.set_yticklabels(y_labels, fontsize=10)
								ax_forest.set_ylim(min(y_ticks) - 0.5, max(y_ticks) + 0.5)
							

					
						# Grid
						ax_forest.grid(axis='x', alpha=0.25, linestyle='--', linewidth=0.6, zorder=0)
						ax_forest.set_axisbelow(True)
					
						# Spines
						ax_forest.spines['top'].set_visible(False)
						ax_forest.spines['right'].set_visible(False)
						ax_forest.spines['left'].set_linewidth(1.2)
						ax_forest.spines['bottom'].set_linewidth(1.2)
					
						# Set white background
						ax_forest.set_facecolor('white')
						ax_forest.set_xlabel('Contribution Magnitude', fontsize=11, weight='bold')
					
						# Add risk/protective labels for first row only
						if idx == 0:
							ax_forest.set_title('Functional Cluster Contributions', fontsize=13, weight='bold', pad=10)
							y_top = max(y_ticks) if y_ticks else 0
							ax_forest.text(xlim[1]*0.95, y_top + 0.8, 'Risk →', 
														ha='right', va='center', fontsize=11,
														weight='bold', style='italic', color='darkred')
							ax_forest.text(xlim[0]*0.95, y_top + 0.8, '← Protective', 
														ha='left', va='center', fontsize=11,
														weight='bold', style='italic', color='darkblue')
							
				# ===== ROW TITLE =====
				# Add subtype label on the left side of the row
				fig.text(0.01, 0.92 - (idx * 0.184), subtype,
								fontsize=18, weight='bold', 
								color='black',  # Changed to black
								va='center', ha='left',
								bbox=dict(boxstyle='round,pad=0.5', 
													facecolor='white', 
													edgecolor='black',
													linewidth=2.5))
			
		# Add legend at the bottom
		legend_elements = [
				plt.Line2D([0], [0], marker=COHORT_MARKERS[cohort], 
									color='w', 
									markerfacecolor=COHORT_COLORS[cohort],
									markeredgecolor='black',
									markeredgewidth=1.0,
									markersize=11, 
									label=f'{cohort}',
									linestyle='None')
				for cohort in cohorts
		]
	
		fig.legend(handles=legend_elements, 
#							loc='best',
							ncol=3,
							frameon=True,
							fancybox=True,
							shadow=True,
							fontsize=12,
							title='ePRS model',
							title_fontsize=13,
							bbox_to_anchor=(0.5, 0.005))
	
		# Save figure
		plt.savefig(output_path, dpi=300, bbox_inches='tight', 
								facecolor='white', edgecolor='none')
		print(f"Combined visualization saved: {output_path}")
		plt.close()
	
	

# ============================================================================
# NEW: CLUSTER SUBTYPE ASSIGNMENT FROM TIERED ENRICHMENT CSV
# ============================================================================

# Maps the source prefix in column names to the corresponding column in cluster_subtype_map.csv
SOURCE_TO_MAP_COLUMN = {
    'kim':        'kim_2023',
    'suzuki':     'suzuki_2024',
    'udler':      'udler_2019',
    'smith':      'smith_2024_ma',
    'yaghootkar': 'yaghootkar_2014',
    'hla':        'hla_cluster',   # custom HLA autoimmune cluster (→ SAID)
}

# ── Directional cluster–feature map ──────────────────────────────────────────
# Expected direction of a clinical feature *within a cluster's phenotypic profile*
# +1 = elevated  (feature value is HIGHER in people carrying this cluster's loci)
# -1 = reduced   (feature value is LOWER)
# Sources: Udler 2019, Suzuki 2024, Yaghootkar 2014, Ahlqvist 2018, Kim 2023
# Keys are lowercase_underscore normalized feature names (partial match applied).
CLUSTER_FEATURE_DIRECTION: dict = {
    'beta_cell_1': {
        'fasting_glucose': +1, 'glucose': +1, 'hba1c': +1,
        '2h_glucose': +1, 'glucose_2h': +1,
        'fasting_insulin': -1,
    },
    'beta_cell_2': {
        'fasting_glucose': +1, 'hba1c': +1,
        '2h_glucose': +1, 'glucose_2h': +1,
        'fasting_insulin': -1,
    },
    'proinsulin': {
        'fasting_glucose': +1, 'hba1c': +1, 'proinsulin': +1,
        'fasting_insulin': -1, 'c_peptide': -1,
    },
    'hyper_insulin': {
        'fasting_insulin': +1, 'c_peptide': +1,
        'triglycerides': +1,
        'hdl': -1, 'hdl_cholesterol': -1,
    },
    'res_glycaemic': {
        'fasting_glucose': +1, 'hba1c': +1,
        # NOT strongly associated with insulin direction
    },
    'obesity': {
        'bmi': +1, 'body_mass_index': +1, 'weight': +1,
        'waist_circumference': +1, 'waist_hip_ratio': +1, 'whr': +1,
        'hip_circumference': +1,
        'body_fat': +1, 'body_fat_percentage': +1,
        'basal_metabolic_rate': +1,
        'hdl': -1, 'hdl_cholesterol': -1,
    },
    'body_fat': {
        # Elevated body fat / visceral fat; NOT strongly tied to BMI or standard lipids
        'body_fat': +1, 'body_fat_percentage': +1,
        'android_fat_mass': +1, 'asat': +1,
        'visceral_adipose_tissue': +1, 'vat': +1,
    },
    'lipodystrophy': {
        # Paradoxical: insulin-resistant despite low fat; HIGH insulin, TG, BP, WHR; LOW fat/HDL
        'fasting_insulin': +1, 'c_peptide': +1,
        'triglycerides': +1,
        'waist_hip_ratio': +1, 'whr': +1,
        'systolic_blood_pressure': +1, 'diastolic_blood_pressure': +1,
        'body_fat': -1, 'body_fat_percentage': -1,
        'bmi': -1, 'body_mass_index': -1,
        'hdl': -1, 'hdl_cholesterol': -1,
        'gfat': -1,   # gluteofemoral adipose tissue
    },
    'metabolic_syndrome': {
        'triglycerides': +1,
        'systolic_blood_pressure': +1, 'diastolic_blood_pressure': +1,
        'fasting_glucose': +1,
        'waist_circumference': +1,
        'hdl': -1, 'hdl_cholesterol': -1,
    },
    'liver_lipid': {
        'alanine_aminotransferase': +1, 'alt': +1,
        'aspartate_aminotransferase': +1, 'ast': +1,
        'gamma_glutamyltransferase': +1, 'ggt': +1,
        'alkaline_phosphatase': +1, 'alp': +1,
        'urate': +1, 'triglycerides': +1,
        # Lower cholesterol transport (lipo A pathway)
        'ldl': -1, 'ldl_cholesterol': -1,
        'total_cholesterol': -1,
    },
}


def _feature_concordance(feature: str, cluster_key: str,
                          effect_size_r, odds_ratio) -> int:
    """
    Compare the observed direction of a feature to the cluster's expected direction.

    Returns
    -------
    +1  concordant  — observed direction matches cluster phenotype
    -1  discordant  — opposite direction
     0  neutral     — no directional information available

    Rules
    -----
    * Genomic loci (have odds_ratio, no effect_size_r): already filtered to OR>1
      (risk direction), so always concordant (+1).
    * Clinical features: compare sign of effect_size_r to CLUSTER_FEATURE_DIRECTION.
    * If the feature is not in the map, return 0 (neutral).
    """
    has_r  = pd.notna(effect_size_r) and float(effect_size_r) != 0.0
    has_or = pd.notna(odds_ratio)

    # Genomic locus: all retained loci are OR>1, always concordant
    if has_or and not has_r:
        return 1

    # Clinical feature without valid r: neutral
    if not has_r:
        return 0

    cluster_map = CLUSTER_FEATURE_DIRECTION.get(cluster_key, {})
    if not cluster_map:
        return 0

    # Normalise feature key: lowercase, underscores
    feat_key = re.sub(r'[\s\-]+', '_', str(feature).lower().strip())

    # Exact match first
    expected = cluster_map.get(feat_key, None)

    # Partial substring match (e.g. 'hdl_direct' matches key 'hdl')
    if expected is None:
        for k, v in cluster_map.items():
            if k in feat_key or feat_key in k:
                expected = v
                break

    if expected is None:
        return 0

    observed = 1 if float(effect_size_r) > 0 else -1
    return 1 if observed == expected else -1


def load_cluster_subtype_map(mapping_file):
    """
    Load the research-based cluster-to-subtype mapping from cluster_subtype_map.csv.

    The file has columns: cluster_key, display_name, subtypes (semicolon-separated),
    udler_2019, kim_2023, smith_2024_ma, yaghootkar_2014, suzuki_2024, assign_manually, ...

    Returns
    -------
    cluster_key_to_subtypes : dict  {cluster_key: [subtype, ...]}
    source_cluster_to_key   : dict  {(source_prefix, cluster_name_lower): cluster_key}
    map_df                  : DataFrame  (full table for downstream use)
    """
    df = pd.read_csv(mapping_file)

    # cluster_key → list[subtype]
    cluster_key_to_subtypes = {}
    for _, row in df.iterrows():
        ck = row['cluster_key']
        subtypes_raw = str(row.get('subtypes', ''))
        subtypes = [s.strip() for s in subtypes_raw.split(';')
                    if s.strip() and s.strip().lower() != 'nan']
        cluster_key_to_subtypes[ck] = subtypes

    # (source_prefix, cluster_name_lower) → cluster_key
    source_cluster_to_key = {}
    for _, row in df.iterrows():
        ck = row['cluster_key']
        for src_prefix, col_name in SOURCE_TO_MAP_COLUMN.items():
            val = row.get(col_name, '')
            if pd.notna(val) and str(val).strip():
                name_lower = str(val).strip().lower()
                source_cluster_to_key[(src_prefix, name_lower)] = ck

    print(f"Loaded cluster_subtype_map: {len(df)} clusters, "
          f"{len(source_cluster_to_key)} source-name entries")
    return cluster_key_to_subtypes, source_cluster_to_key, df


def map_tiered_columns_to_cluster_keys(tiered_df, source_cluster_to_key):
    """
    Map each binary cluster column in the tiered enrichment CSV to a cluster_key.

    Column format: {source}_functional_cluster_{loci|clinical}_map__{ClusterName}
    Source prefix: first token before the first underscore (kim, suzuki, udler, smith).

    Returns
    -------
    col_to_cluster_key : dict  {column_name: cluster_key or None}
    """
    cluster_cols = [
        c for c in tiered_df.columns
        if '__' in c and any(c.startswith(p + '_') for p in SOURCE_TO_MAP_COLUMN)
    ]

    import re

    def _normalise(s):
        """Lower-case; replace hyphens/slashes/underscores with space;
        remove filler word 'and'; insert space before digits after letters
        ('cell1' → 'cell 1'); collapse multiple spaces."""
        s = s.lower()
        s = re.sub(r'[-/_]', ' ', s)
        s = re.sub(r'\band\b', ' ', s)            # 'liver and lipid' → 'liver  lipid'
        s = re.sub(r'([a-z])(\d)', r'\1 \2', s)   # 'cell1' → 'cell 1'
        s = re.sub(r'\s+', ' ', s).strip()
        return s

    # Pre-normalise the lookup keys
    normalised_lookup = {
        (src, _normalise(name)): ck
        for (src, name), ck in source_cluster_to_key.items()
    }

    col_to_cluster_key = {}
    unmapped = []

    for col in cluster_cols:
        # Skip pandas unnamed artifact columns
        if 'Unnamed' in col:
            col_to_cluster_key[col] = None
            continue

        source_prefix = col.split('_')[0]
        raw_name = col.split('__', 1)[1].strip()
        cluster_name = _normalise(raw_name)

        # Exact match on normalised string
        key = normalised_lookup.get((source_prefix, cluster_name))

        # Fallback: partial match on normalised strings
        if key is None:
            for (src, norm_name), ck in normalised_lookup.items():
                if src == source_prefix and (norm_name in cluster_name or cluster_name in norm_name):
                    key = ck
                    break

        col_to_cluster_key[col] = key
        if key is None:
            unmapped.append(col)

    if unmapped:
        print(f"  Warning: {len(unmapped)} cluster columns could not be mapped to a cluster_key:")
        for col in unmapped:
            print(f"    {col}")

    mapped_count = sum(1 for v in col_to_cluster_key.values() if v is not None)
    print(f"  Mapped {mapped_count}/{len(cluster_cols)} cluster columns to cluster_keys")
    return col_to_cluster_key


def assign_subtypes_to_features(tiered_df, col_to_cluster_key,
                                 cluster_key_to_subtypes, cluster_subtype_map_df,
                                 cohort_col='cohort'):
    """
    For each feature row in the tiered enrichment CSV, assign subtypes based on
    which functional cluster binary columns equal 1.

    A feature in cohort C mapped to a cluster that links to subtype S produces one
    record (feature, C, cluster_key, S).  If the same gene appears in two cohorts
    it generates two records, one per cohort — as requested.

    Returns
    -------
    assigned_df : long-format DataFrame with one row per feature×cohort×cluster_key×subtype
    """
    key_to_display  = dict(zip(cluster_subtype_map_df['cluster_key'],
                                cluster_subtype_map_df['display_name']))
    key_to_manual   = dict(zip(cluster_subtype_map_df['cluster_key'],
                                cluster_subtype_map_df['assign_manually'].fillna(0).astype(int)))

    valid_cluster_cols = {col: ck for col, ck in col_to_cluster_key.items()
                          if ck is not None}

    meta_cols = [
        'feature', 'cohort', 'training_data_used', 'feature_source',
        'odds_ratio', 'effect_size_r', 'specificity_tier',
        'n_exclusive_cohorts', 'exclusive_cohorts',
        'snp1_gene_names', 'snp2_gene_names', 'mahajan_clinical_feature',
    ]
    present_meta = [c for c in meta_cols if c in tiered_df.columns]

    rows = []
    for _, row in tiered_df.iterrows():
        cohort = row.get(cohort_col)
        if pd.isna(cohort):
            continue

        for col, cluster_key in valid_cluster_cols.items():
            val = row.get(col)
            if pd.notna(val) and val == 1:
                conc = _feature_concordance(
                    row.get('feature', ''),
                    cluster_key,
                    row.get('effect_size_r'),
                    row.get('odds_ratio'),
                )
                for subtype in cluster_key_to_subtypes.get(cluster_key, []):
                    entry = {c: row.get(c) for c in present_meta}
                    entry['cluster_key']    = cluster_key
                    entry['display_name']   = key_to_display.get(cluster_key, cluster_key)
                    entry['subtype']        = subtype
                    entry['assign_manually'] = key_to_manual.get(cluster_key, 0)
                    entry['source_column']  = col
                    entry['concordance']    = conc
                    rows.append(entry)

    assigned_df = pd.DataFrame(rows)
    print(f"  Assigned {len(assigned_df)} feature×cohort×cluster×subtype records "
          f"({assigned_df['feature'].nunique() if len(assigned_df) else 0} unique features)")
    return assigned_df


def compute_subtype_cohort_enrichment(assigned_df, tiered_df,
                                       col_to_cluster_key, cohort_col='cohort'):
    """
    Compute per-subtype cohort enrichment with two normalizations.

    Normalization 1 — by cohort feature total
        norm_by_cohort = n_features(subtype, cohort) / total_features(cohort)

    Normalization 2 — by cluster gene count per cohort
        norm_by_cluster = n_features(subtype, cluster, cohort)
                          / total_features_in_cluster(cluster, cohort)

    Returns
    -------
    subtype_cohort_df   : DataFrame  enrichment summary (subtype × cohort)
    subtype_cluster_df  : DataFrame  detailed breakdown (subtype × cluster × cohort)
    """
    # --- Total features per cohort in the tiered CSV ---
    cohort_totals = tiered_df[cohort_col].dropna().value_counts().to_dict()

    # --- Features per cluster_key per cohort in the tiered CSV ---
    # De-duplicate by (feature, cohort, cluster_key) so a feature matched by both
    # clinical and loci map for the same cluster only counts once
    valid_cols = {col: ck for col, ck in col_to_cluster_key.items()
                  if ck is not None and col in tiered_df.columns}

    cluster_cohort_totals = {}   # {cluster_key: {cohort: count}}
    seen = set()
    for col, cluster_key in valid_cols.items():
        sub = tiered_df[tiered_df[col] == 1][[cohort_col, 'feature']]
        for _, r in sub.iterrows():
            cohort = r[cohort_col]
            feature = r['feature']
            if pd.isna(cohort):
                continue
            dedup_key = (cluster_key, cohort, feature)
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            cluster_cohort_totals.setdefault(cluster_key, {})
            cluster_cohort_totals[cluster_key][cohort] = \
                cluster_cohort_totals[cluster_key].get(cohort, 0) + 1

    # --- De-duplicate assigned_df before counting ---
    dedup = assigned_df.drop_duplicates(
        subset=['feature', 'cohort', 'cluster_key', 'subtype']
    )

    # --- Normalization 1: subtype × cohort ---
    sc = (dedup.groupby(['subtype', 'cohort']).size().reset_index(name='n_features'))
    sc['total_cohort_features'] = sc['cohort'].map(cohort_totals).fillna(0)
    sc['norm_by_cohort_features'] = (
        sc['n_features'] / sc['total_cohort_features'].replace(0, np.nan)
    )

    # --- Normalization 2: subtype × cluster × cohort ---
    scc = (dedup.groupby(['subtype', 'cluster_key', 'cohort'])
               .size().reset_index(name='n_features_in_cluster'))

    def _cluster_cohort_total(row):
        return cluster_cohort_totals.get(row['cluster_key'], {}).get(row['cohort'], 0)

    scc['n_cluster_genes_in_cohort'] = scc.apply(_cluster_cohort_total, axis=1)
    scc['norm_by_cluster_gene_count'] = (
        scc['n_features_in_cluster'] /
        scc['n_cluster_genes_in_cohort'].replace(0, np.nan)
    )

    return sc, scc


# ============================================================================
# DIRECT EVIDENCE: subtype*_map.csv FILES FROM LITERATURE
# ============================================================================

# Recognised subtype names (used to detect wide-format files)
KNOWN_SUBTYPES = {'SAID', 'SIDD', 'SIRD', 'MOD', 'MARD', 'LADA', 'MODY'}


def load_subtype_direct_maps(subtype_research_path: str) -> List[dict]:
    """
    Scan *subtype_research_path* for files matching ``subtype*_map.csv``
    (case-insensitive glob) and parse each into feature→subtype assignments.

    Two formats are handled automatically:

    **Wide format** (feature as first / unlabelled column, subtypes as other columns)::

        feature    SAID  SIDD  SIRD  MOD  MARD
        BMI        23.1  24.7  28.3  26.1  24.0
        HbA1c      70.5  89.1  57.4  53.0  54.1

    For each feature, the subtype(s) whose column value has a Z-score ≥ *z_thresh*
    (default 0.5) above the cross-subtype mean are considered associated.
    At minimum the top-scoring subtype is always returned.

    **Long format** (one row per feature×subtype)::

        feature     subtype  effect_size
        rs7903146   SIDD     0.8
        rs7903146   SAID     0.4

    Returns
    -------
    list of dicts, each with keys:
        source_file  : str   basename of the CSV
        feature      : str   normalised feature name
        subtype      : str   T2D subtype
        evidence     : float column value (or 1.0 for long-format presence)
        z_score      : float Z-score across subtypes (NaN for long format)
    """
    import glob as _glob

    pattern = os.path.join(subtype_research_path, 'subtype*_map.csv')
    files = sorted(_glob.glob(pattern, recursive=False))
    if not files:
        return []

    print(f"\n  Found {len(files)} subtype direct-map file(s):")
    all_records = []

    for fpath in files:
        src = os.path.basename(fpath)
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            print(f"    [WARN] Could not load {src}: {e}")
            continue

        # Detect subtype columns — any column whose name (uppercased) is in KNOWN_SUBTYPES
        subtype_cols = [c for c in df.columns if c.strip().upper() in KNOWN_SUBTYPES]

        if subtype_cols:
            # ---- Wide format ----
            # Feature column = first column not in subtype_cols
            feat_col = next((c for c in df.columns if c not in subtype_cols), None)
            if feat_col is None:
                print(f"    [WARN] {src}: no feature column found — skipping")
                continue

            n_recs = 0
            for _, row in df.iterrows():
                feat = str(row[feat_col]).strip()
                if not feat or feat.lower() == 'nan':
                    continue

                vals = pd.to_numeric(
                    pd.Series({c: row[c] for c in subtype_cols}), errors='coerce'
                ).dropna()
                if vals.empty:
                    continue

                # Z-score across subtypes for this feature
                mean_v, std_v = vals.mean(), vals.std()
                z_scores = (vals - mean_v) / std_v if std_v > 0 else pd.Series(0.0, index=vals.index)

                z_thresh = 0.5
                associated = z_scores[z_scores >= z_thresh].index.tolist()
                # Always include the maximum subtype even below threshold
                if not associated:
                    associated = [z_scores.idxmax()]

                for sub in associated:
                    all_records.append({
                        'source_file': src,
                        'feature':     feat,
                        'subtype':     sub.strip().upper(),
                        'evidence':    float(vals[sub]),
                        'z_score':     float(z_scores[sub]),
                    })
                    n_recs += 1

            print(f"    {src}  [wide]: {len(df)} features × {len(subtype_cols)} subtypes "
                  f"→ {n_recs} feature×subtype records")

        elif 'subtype' in df.columns:
            # ---- Long format ----
            feat_candidates = [c for c in df.columns
                               if c.lower() in ('feature', 'locus', 'loci', 'snp', 'rsid',
                                                 'gene', 'cluster', 'cluster_key')]
            feat_col = feat_candidates[0] if feat_candidates else df.columns[0]
            ev_col   = next((c for c in df.columns
                             if c.lower() in ('effect_size', 'beta', 'or', 'odds_ratio',
                                              'z_score', 'score', 'value')), None)

            n_recs = 0
            for _, row in df.iterrows():
                feat = str(row[feat_col]).strip()
                sub  = str(row['subtype']).strip().upper()
                if not feat or feat.lower() == 'nan' or sub not in KNOWN_SUBTYPES:
                    continue
                ev = float(row[ev_col]) if ev_col and pd.notna(row.get(ev_col)) else 1.0
                all_records.append({
                    'source_file': src,
                    'feature':     feat,
                    'subtype':     sub,
                    'evidence':    ev,
                    'z_score':     float('nan'),
                })
                n_recs += 1

            print(f"    {src}  [long]: {n_recs} feature×subtype records")

        else:
            print(f"    [WARN] {src}: no recognised subtype columns or 'subtype' column — skipping")

    return all_records


def match_direct_maps_to_tiered(direct_records: List[dict],
                                 tiered_df: pd.DataFrame,
                                 cohort_col: str = 'cohort') -> pd.DataFrame:
    """
    Match features from ``direct_records`` (output of
    :func:`load_subtype_direct_maps`) against rows in *tiered_df*.

    Matching priority (first hit wins, case-insensitive):
    1. ``feature`` column in tiered_df  (exact, normalised)
    2. ``mahajan_clinical_feature`` column  (clinical bridge)
    3. ``snp1_gene_names`` / ``snp2_gene_names``  (gene-symbol substring)

    Each matched tiered row produces one record per (matched_feature × cohort × subtype).

    Returns
    -------
    DataFrame compatible with ``assigned_df`` from
    :func:`assign_subtypes_to_features`.
    """
    if not direct_records:
        return pd.DataFrame()

    def _norm(s):
        return str(s).strip().lower() if pd.notna(s) else ''

    # Build lookup indices on tiered_df
    feat_idx     = {_norm(v): idxs
                    for v, idxs in tiered_df.groupby(tiered_df['feature'].map(_norm)).groups.items()}
    mahajan_col  = 'mahajan_clinical_feature'
    mahajan_idx  = {}
    if mahajan_col in tiered_df.columns:
        for v, idxs in tiered_df.groupby(tiered_df[mahajan_col].map(_norm)).groups.items():
            if v:
                mahajan_idx[v] = idxs

    rows = []
    for rec in direct_records:
        raw = rec['feature']
        norm_feat = _norm(raw)

        # Try exact feature match first
        matched_indices = feat_idx.get(norm_feat, [])

        # Try Mahajan clinical bridge
        if len(matched_indices) == 0 and mahajan_idx:
            matched_indices = mahajan_idx.get(norm_feat, [])

        # HLA gene-region prefix match:
        # Map entries like 'HLA-DQB1' match tiered features 'HLA-dqb1*504', 'HLA-dqb1*301' etc.
        # Pattern: map feature = 'HLA-<GENE>' (no allele suffix), tiered feature = 'HLA-<gene>*<allele>'
        if len(matched_indices) == 0 and re.match(r'^hla-[a-z0-9]+$', norm_feat):
            matched_indices = tiered_df[
                tiered_df['feature'].str.lower().str.startswith(norm_feat, na=False)
            ].index.tolist()

        # Try gene name substring (token match in snp1_gene_names / snp2_gene_names).
        # Only applied to features that look like genomic identifiers (gene symbols,
        # SNP IDs, locus names) — NOT to clinical trait names like "eGFR", "BMI",
        # "Glucose" which could accidentally match gene abbreviations.
        def _looks_genomic(s):
            """Return True if s resembles a gene symbol, SNP ID, or locus name."""
            s = str(s).strip()
            if re.match(r'^rs\d+$', s, re.IGNORECASE):          # rsXXXX SNP ID
                return True
            if re.match(r'^HLA[-_*]', s, re.IGNORECASE):        # HLA alleles/genes
                return True
            # Gene symbols: all uppercase + digits/hyphens, ≥ 4 chars, no spaces
            if re.match(r'^[A-Z][A-Z0-9_-]{3,}$', s):
                return True
            return False

        if len(matched_indices) == 0 and _looks_genomic(raw):
            gene_tokens = set(re.split(r'[;,\s]+', norm_feat)) - {''}
            for g in gene_tokens:
                if not g:
                    continue
                for col in ('snp1_gene_names', 'snp2_gene_names'):
                    if col not in tiered_df.columns:
                        continue
                    mask = tiered_df[col].str.contains(g, case=False, na=False, regex=False)
                    matched_indices = tiered_df[mask].index.tolist()
                    if matched_indices:
                        break
                if matched_indices:
                    break

        for idx in matched_indices:
            trow = tiered_df.loc[idx]
            cohort = trow.get(cohort_col)
            if pd.isna(cohort):
                continue
            rows.append({
                'feature':         trow.get('feature', raw),
                'cohort':          cohort,
                'subtype':         rec['subtype'],
                'cluster_key':     f"direct__{rec['source_file'].replace('.csv', '')}",
                'display_name':    f"Direct ({rec['source_file']})",
                'assign_manually': 0,
                'source_column':   rec['source_file'],
                'training_data_used': trow.get('training_data_used'),
                'feature_source':  trow.get('feature_source'),
                'odds_ratio':      trow.get('odds_ratio'),
                'effect_size_r':   trow.get('effect_size_r'),
                'specificity_tier': trow.get('specificity_tier'),
                'evidence':        rec['evidence'],
                'z_score_direct':  rec['z_score'],
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        print(f"\n  Direct-map matches: {len(result)} feature×cohort×subtype records "
              f"({result['feature'].nunique()} unique features)")
    return result


if __name__ == "__main__":

    # ------------------------------------------------------------------
    # CLI argument parsing
    # ------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description=(
            "Cluster-subtype assignment from tiered enrichment CSV. "
            "Joins cluster_subtype_map.csv onto functional_cluster_enrichment_tiered.csv, "
            "assigns subtypes to each feature×cohort, and computes normalised enrichment."
        )
    )
    parser.add_argument(
        '--tiered_csv', required=True,
        help='Path to functional_cluster_enrichment_tiered.csv'
    )
    parser.add_argument(
        '--cluster_subtype_map', required=True,
        help='Path to cluster_subtype_map.csv (subtypeResearch/)'
    )
    parser.add_argument(
        '--out_dir', required=True,
        help='Output directory for the annotated and enrichment CSVs'
    )
    parser.add_argument(
        '--cluster_weights', default=None,
        help='(Optional) Path to clusterWeights.csv for the legacy OR-weighted analysis'
    )
    parser.add_argument(
        '--legacy_mapping', default=None,
        help='(Optional) Path to clusterToSubtype.csv for the legacy analysis'
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 80)
    print("SUBTYPE ASSIGNMENT FROM TIERED ENRICHMENT CSV")
    print("=" * 80)

    # ------------------------------------------------------------------
    # Step 1: Load tiered enrichment CSV
    # ------------------------------------------------------------------
    print(f"\nLoading tiered enrichment CSV: {args.tiered_csv}")
    tiered_df = pd.read_csv(args.tiered_csv)
    print(f"  {len(tiered_df)} rows, {len(tiered_df.columns)} columns")
    print(f"  Cohorts: {sorted(tiered_df['cohort'].dropna().unique().tolist())}")

    # ------------------------------------------------------------------
    # Step 2: Load cluster-subtype map
    # ------------------------------------------------------------------
    print(f"\nLoading cluster-subtype map: {args.cluster_subtype_map}")
    cluster_key_to_subtypes, source_cluster_to_key, cluster_map_df = \
        load_cluster_subtype_map(args.cluster_subtype_map)

    # ------------------------------------------------------------------
    # Step 3: Map tiered CSV columns → cluster_keys
    # ------------------------------------------------------------------
    print("\nMapping cluster columns to cluster_keys...")
    col_to_cluster_key = map_tiered_columns_to_cluster_keys(tiered_df, source_cluster_to_key)

    # ------------------------------------------------------------------
    # Step 4: Assign subtypes to features
    # ------------------------------------------------------------------
    print("\nAssigning subtypes to features...")
    assigned_df = assign_subtypes_to_features(
        tiered_df, col_to_cluster_key,
        cluster_key_to_subtypes, cluster_map_df
    )

    # ------------------------------------------------------------------
    # Step 4b: Direct-evidence from subtype*_map.csv files
    # ------------------------------------------------------------------
    subtype_research_path = os.path.dirname(args.cluster_subtype_map)
    print(f"\nScanning for subtype direct-map files in: {subtype_research_path}")
    direct_records = load_subtype_direct_maps(subtype_research_path)

    if direct_records:
        direct_df = match_direct_maps_to_tiered(direct_records, tiered_df)
        if not direct_df.empty:
            assigned_df = pd.concat([assigned_df, direct_df], ignore_index=True)
            print(f"  Combined assigned_df: {len(assigned_df)} total records")
            # Save direct-evidence breakdown separately for inspection
            direct_path = os.path.join(args.out_dir, 'subtype_direct_evidence.csv')
            direct_df.to_csv(direct_path, index=False)
            print(f"  Saved direct evidence → {direct_path}")
    else:
        print("  No subtype*_map.csv files found — skipping direct-evidence step.")

    # Save long-format annotated table
    annotated_path = os.path.join(args.out_dir, 'subtype_feature_assignments.csv')
    assigned_df.to_csv(annotated_path, index=False)
    print(f"\n  Saved annotated assignments → {annotated_path}")

    # ------------------------------------------------------------------
    # Step 5: Compute enrichment (two normalizations)
    # ------------------------------------------------------------------
    print("\nComputing subtype-cohort enrichment...")
    subtype_cohort_df, subtype_cluster_df = compute_subtype_cohort_enrichment(
        assigned_df, tiered_df, col_to_cluster_key
    )

    # Save enrichment tables
    enrich_cohort_path = os.path.join(args.out_dir, 'subtype_cohort_enrichment.csv')
    enrich_cluster_path = os.path.join(args.out_dir, 'subtype_cluster_cohort_enrichment.csv')
    subtype_cohort_df.to_csv(enrich_cohort_path, index=False)
    subtype_cluster_df.to_csv(enrich_cluster_path, index=False)
    print(f"  Saved subtype×cohort enrichment     → {enrich_cohort_path}")
    print(f"  Saved subtype×cluster×cohort detail → {enrich_cluster_path}")

    # ------------------------------------------------------------------
    # Step 6: Console summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("SUBTYPE COHORT ENRICHMENT SUMMARY")
    print("=" * 80)
    print("\nn_features per subtype×cohort:")
    pivot_n = subtype_cohort_df.pivot(
        index='subtype', columns='cohort', values='n_features'
    ).fillna(0).astype(int)
    print(pivot_n.to_string())

    print("\nnorm_by_cohort_features (fraction of cohort's features assigned to subtype):")
    pivot_norm = subtype_cohort_df.pivot(
        index='subtype', columns='cohort', values='norm_by_cohort_features'
    ).fillna(0)
    print(pivot_norm.round(4).to_string())

    if len(assigned_df):
        print("\nassign_manually=1 entries (review recommended):")
        manual = assigned_df[assigned_df['assign_manually'] == 1][
            ['feature', 'cohort', 'cluster_key', 'subtype']
        ].drop_duplicates()
        if len(manual):
            print(manual.to_string(index=False))
        else:
            print("  None — all assignments are automatic.")

    # ------------------------------------------------------------------
    # Step 7: (Optional) legacy OR-weighted analysis
    # ------------------------------------------------------------------
    if args.cluster_weights and args.legacy_mapping:
        pheno_path = os.path.dirname(os.path.dirname(args.cluster_weights))
        print("\n" + "=" * 80)
        print("LEGACY OR-WEIGHTED ANALYSIS (clusterWeights.csv)")
        print("=" * 80)

        df = pd.read_csv(args.cluster_weights)
        results = rank_functional_clusters_by_cohort(df)
        display_direction_analysis(results)

        cluster_mapping = load_cluster_subtype_mapping(args.legacy_mapping)

        subtype_weighted = calculate_weighted_subtype_enrichment(
            df, results, cluster_mapping, weighting='effect_size'
        )
        display_weighted_results(subtype_weighted)

        fig_dir = os.path.join(pheno_path, 'figures', 'clinicalFigures')
        os.makedirs(fig_dir, exist_ok=True)

        forest_data = prepare_forest_plot_data(
            subtype_weighted['risk_contributions'],
            subtype_weighted['protective_contributions'],
            results['risk_scores'],
            cluster_mapping
        )
        create_combined_visualization(
            subtype_weighted['cohort_percentages_per_subtype'],
            forest_data,
            os.path.join(fig_dir, 'combined_pie_forest_weighted_by_effect_size.riskOnly.png')
        )
        plot_weighted_subtype_pies(
            subtype_weighted,
            os.path.join(fig_dir, 'pie_weighted_by_effect_size.png')
        )
        plot_risk_protective_comparison(
            results,
            save_path=os.path.join(pheno_path, 'figures', 'risk_protective_comparison.png')
        )

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print("\nOutputs written to:", args.out_dir)