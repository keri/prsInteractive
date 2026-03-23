#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import re




def compare_gene_env_to_genetic(inputFile, config_file):
    df = pd.read_csv(inputFile)
    
    cardio_features = df[
        (df['model'] == 'cardio') & (df['feature'].str.contains(','))
    ]['feature'].tolist()
    
    mainDf = df[df['model'] == 'all+env_combined']
    
    # ── Bug 1 fix: start dfFiltered from the original df and update it
    #    incrementally so each model's filtering is preserved.
    dfFiltered = df.copy()
    
    for model in ['cardio', 'all']:
        
        gene_env = df[
            (df['model'] == model) &
            (df['coefs'] != 0) &
            (df['feature'].isin(cardio_features))
        ]
        n_input = len(gene_env)
        
        filterList = []
        kept_more_extreme = []
        
        print(f"\nAnalyzing {n_input} GxGxE features...")
        
        for feature in gene_env['feature']:
            g = ','.join(feature.split(',')[1:])
            e = feature.split(',')[0]
            
            eDf = mainDf[mainDf['feature'] == feature]
            eValue = eDf['coefs'].values[0] if not eDf.empty else 0
            
            gDf = mainDf[mainDf['feature'] == g]
            mDf = mainDf[mainDf['feature'] == e]
            
            if gDf.empty:
                print(f'no corresponding genetic feature found : {g}')
                gValue = 0
            else:
                gValue = gDf['coefs'].values[0]
                
            if mDf.empty:
                print(f'environmental feature {e} not found, setting to 0')
                mValue = 0
            else:
                mValue = mDf['coefs'].values[0]
                
            same_direction = (
                (eValue > 0 and gValue >= 0 and mValue >= 0) or
                (eValue < 0 and gValue <= 0 and mValue <= 0)
            )
            
            if same_direction:
                # All three in the same direction: use directional comparison
                if eValue > 0:
                    keep = (eValue >= gValue) and (eValue >= mValue)
                else:
                    keep = (eValue <= gValue) and (eValue <= mValue)
            else:
                # Opposite directions: keep only if the interaction is more
                # extreme in magnitude than BOTH components individually
                keep = (abs(eValue) > abs(gValue)) and (abs(eValue) > abs(mValue))
                
            if keep:
                kept_more_extreme.append(feature)
                print(f'kept GxGxE feature :  {feature}')
                print(f'which is more extreme than {g} : {gValue}')
                print(f'and more extreme than {e} : {mValue}')
            else:
                filterList.append(feature)
                
        #    so the previous model's filtering is not lost.
        mask = (dfFiltered['model'] == model) & (dfFiltered['feature'].isin(filterList))
        dfFiltered = dfFiltered[~mask]
        
        #    consistent with how the input count is defined.
        n_output = len(
            dfFiltered[
                (dfFiltered['model'] == model) &
                (dfFiltered['coefs'] != 0) &
                (dfFiltered['feature'].isin(cardio_features))
            ]
        )
        
        print(f'\n{"="*70}')
        print('Filtered non-additive GxGxE features:')
        print(f'  Input file had:  {n_input} GxGxE features')
        print(f'  Filtered out:    {len(filterList)} features (not more extreme than GxG)')
        print(f'  Kept extreme:    {len(kept_more_extreme)} features (more extreme than GxG)')
        print(f'  Output file has: {n_output} GxGxE features')
        print(f'{"="*70}\n')
        
    # Save once, after both models have been processed
    inputFilePrefix = inputFile.rsplit('.', 1)[0]
    filteredFile = f'{inputFilePrefix}.filtered.csv'
    dfFiltered.to_csv(filteredFile, index=False)
    print(f'  Filtered file saved to: {filteredFile}')
    
    # Update config file (unchanged)
    try:
        with open(config_file, 'r') as f:
            lines = f.readlines()
        updated_lines = [
            l for l in lines
            if not l.strip().startswith('FEATURE_SCORES_FILE_FILTERED=')
        ]
        updated_lines.append(f'FEATURE_SCORES_FILE_FILTERED="{filteredFile}"\n')
        with open(config_file, 'w') as f:
            f.writelines(updated_lines)
        print(f'Config file updated: {config_file}')
    except FileNotFoundError:
        print(f'Warning: Config file not found: {config_file}')
    except Exception as ex:
        print(f'Warning: Could not update config file: {ex}')
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="filtering non-additive features in cardio model ....")
    parser.add_argument("--feature_scores_file", help="Path to feature scores file post modelling folder")
    parser.add_argument("--config_file", help="Path to pheno.config")
    
    args = parser.parse_args()

    
    # Prefer command-line input if provided; fallback to env var
    feature_scores_file = args.feature_scores_file or os.environ.get('FEATURE_SCORES_FILE')
    config_file = args.config_file or os.environ.get('CONFIG_FILE')
    
#   pheno = 'type2Diabetes'
#   feature_scores_file = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi/scores/featureScoresReducedFinalModel.csv'
#   config_file = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi/pheno.config'
        
    compare_gene_env_to_genetic(feature_scores_file,config_file)