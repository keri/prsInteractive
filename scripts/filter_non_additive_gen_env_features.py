#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import re



def compare_gene_env_to_genetic(inputFile, config_file):
    """
    Filter GxGxE features to keep only those that are MORE extreme than their GxG component.
    
    Logic:
    - For positive coefs (risk-increasing): Keep GxGxE if it's MORE positive than GxG
    - For negative coefs (protective): Keep GxGxE if it's MORE negative than GxG
    
    Parameters:
    -----------
    inputFile : str
        Path to input feature scores CSV
    config_file : str
        Path to config file to update
    """
    df = pd.read_csv(inputFile)
    
    # Count original GxGxE features
    n_cardio = len(df[(df['model'] == 'cardio') & (df['coefs'] != 0) & (df['feature'].str.contains(','))])
    
    # Get non-zero GxGxE features (cardio model with interactions)
    gene_env = df[(df['model'] == 'cardio') & (df['coefs'] != 0) & (df['feature'].str.contains(','))]
    
    filterList = []
    kept_more_extreme = []
    
    print(f"\nAnalyzing {len(gene_env)} GxGxE features...")
    
    for feature in gene_env['feature']:
        # Extract genetic component (everything after first element)
        # E.g., "BMI,rs123,rs456" -> "rs123,rs456"
        g = ','.join(feature.split(',')[1:])
        
        # Get GxGxE coefficient
        eDf = gene_env[gene_env['feature'] == feature]
        eValue = eDf['coefs'].values[0]
        
        # Find corresponding GxG feature in epi or main models
        gDf = df[(df['feature'] == g) & (df['model'].isin(['epi', 'main']))]
        
        if gDf.empty:
            # No corresponding GxG feature found - keep GxGxE
            continue
        else:

            # Filter logic: Keep only if GxGxE is MORE extreme than GxG
            if eValue > 0:  # GxGxE is risk-increasing
                # Get the strongest GxG coefficient if multiple models have it
                gValue = gDf.sort_values(['coefs'], ascending=False).head(1)['coefs'].values[0]
                if gValue > 0: #ensure the coefs are in the same direction
                    # Keep only if GxGxE is MORE positive (more risky)
                    if eValue <= gValue:
                        filterList.append(feature)
                        # print(f"  Filter: {feature[:50]}... | GxG={gValue:.4f}, GxGxE={eValue:.4f} (less risky)")
                    else:
                        kept_more_extreme.append(feature)
                else:
                    filterList.append(feature)
                    
            else:  # GxG is protective (negative)
                # Get the strongest GxG coefficient if multiple models have it
                gValue = gDf.sort_values(['coefs'], ascending=True).head(1)['coefs'].values[0]
                if gValue < 0: #ensure the coefs are in the same direction
                    # Keep only if GxGxE is MORE negative (more protective)
                    if eValue >= gValue:
                        filterList.append(feature)
                        # print(f"  Filter: {feature[:50]}... | GxG={gValue:.4f}, GxGxE={eValue:.4f} (less protective)")
                    else:
                        kept_more_extreme.append(feature)
                else:
                    filterList.append(feature)
                    
    # Filter out non-additive features
    dfFiltered = df[~df['feature'].isin(filterList)]
    
    # Count filtered GxGxE features
    n_filtered = len(dfFiltered[(dfFiltered['model'] == 'cardio') & 
                                 (dfFiltered['coefs'] != 0) & 
                                 (dfFiltered['feature'].str.contains(','))])
    
    # Save filtered file
    inputFilePrefix = inputFile.rsplit('.', 1)[0]  # Handle multiple dots in filename
    filteredFile = f'{inputFilePrefix}.filtered.csv'
    dfFiltered.to_csv(filteredFile, index=False)
    
    # Print summary
    print(f'\n{"="*70}')
    print('Filtered non-additive GxGxE features:')
    print(f'  Input file had:  {n_cardio} GxGxE features')
    print(f'  Filtered out:    {len(filterList)} features (not more extreme than GxG)')
    print(f'  Kept extreme:    {len(kept_more_extreme)} features (more extreme than GxG)')
    print(f'  Output file has: {n_filtered} GxGxE features')
    print(f'  Filtered file saved to: {filteredFile}')
    print(f'{"="*70}\n')
    
    # Update config file
    try:
        with open(config_file, 'r') as f:
            lines = f.readlines()
            
        # Update the FEATURE_SCORES_FILE line
        updated_lines = []
        config_updated = False
        
        for line in lines:
            if line.strip().startswith('FEATURE_SCORES_FILE='):
                updated_lines.append(f'FEATURE_SCORES_FILE="{filteredFile}"\n')
                config_updated = True
            else:
                updated_lines.append(line)
                
        # Write back
        with open(config_file, 'w') as f:
            f.writelines(updated_lines)
            
        if config_updated:
            print(f'Config file updated: {config_file}')
        else:
            print(f'Warning: FEATURE_SCORES_FILE not found in config file')
            
    except FileNotFoundError:
        print(f'Warning: Config file not found: {config_file}')
    except Exception as e:
        print(f'Warning: Could not update config file: {e}')
        
    return dfFiltered, filterList, kept_more_extreme
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="filtering non-additive features in cardio model ....")
    parser.add_argument("--feature_scores_file", help="Path to feature scores file post modelling folder")
    parser.add_argument("--config_file", help="Path to pheno.config")
    
    args = parser.parse_args()

    
    # Prefer command-line input if provided; fallback to env var
    feature_scores_file = args.feature_scores_file or os.environ.get('FEATURE_SCORES_FILE')
    config_file = args.config_file or os.environ.get('CONFIG_FILE')
    
#   feature_scores_file = '/Users/kerimulterer/prsInteractive/results/celiacDisease/summedEpi/scores/featureScoresReducedFinalModel.csv'
#   config_file = '/Users/kerimulterer/prsInteractive/results/celiacDisease/summedEpi/pheno.config'
        
    compare_gene_env_to_genetic(feature_scores_file,config_file)