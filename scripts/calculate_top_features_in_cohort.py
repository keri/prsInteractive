#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
#import pingouin as pg

#import scipy
#import os
#from scipy import stats
#from sklearn.model_selection import cross_val_score, StratifiedKFold
#from sklearn.feature_selection import RFECV
from sklearn.metrics import roc_auc_score, make_scorer, roc_curve, auc
from sklearn.impute import SimpleImputer
#from imblearn.over_sampling import SMOTE
import warnings

from sklearn.preprocessing import StandardScaler

from sklearn.ensemble import RandomForestClassifier

from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import os

from helper.download import *
from helper.data_wrangling import *
from helper.draw_plots import *
from helper.calculate_shap_values import *
#from helper.find_important_features_from_shap import *

# Suppress the specific warning
warnings.filterwarnings("ignore", message="Bad circle positioning")


def find_opposite_sign_features(df, feature_sets,threshold):
    """
    Find features where the matching cohort has positive SHAP z-score 
    and ALL other cohorts have negative z-scores.
    
    Parameters:
    -----------
    df : DataFrame
        Must contain columns: 'features_used', 'cohort', 'feature', 'z_score', 'coef'
    feature_sets : list
        List of feature sets to analyze
    
    Returns:
    --------
    DataFrame with columns: feature_set, feature, matching_z_score, 
                           n_other_cohorts, coef, odds_ratio
    """
    rows = []
    
    for f in feature_sets:
        # Filter for current feature set
        df_feature_set = df[df['features_used'] == f]
        
        # Separate matching cohort from others
        df_matching = df_feature_set[(df_feature_set['cohort'] == f) & (df_feature_set['z_score'] > threshold)]
        df_other = df_feature_set[(df_feature_set['cohort'] != f) & (df_feature_set['feature'].isin(df_matching.feature.tolist()))]
        
        # Get unique features for this matching cohort
        for feature in df_matching['feature'].unique():
            # Get z-score and coef for this feature in matching cohort
            matching_row = df_matching[df_matching['feature'] == feature]
            
            if len(matching_row) == 0:
                continue
            
            matching_z = matching_row['z_score'].values[0]
            matching_coef = matching_row['coefs'].values[0]
            
            # Get z-scores for this feature in all other cohorts
            other_z = df_other[df_other['feature'] == feature]['z_score'].values
            
            # Skip if this feature doesn't appear in other cohorts
            if len(other_z) == 0:
                continue
            
            # Only keep if matching cohort has positive z-score 
            # AND ALL other cohorts have negative z-scores
            if matching_z > 0 and all(z < 0 for z in other_z):
                rows.append({
                    'feature_set': f,
                    'feature': feature,
                    'matching_z_score': matching_z,
                    'n_other_cohorts': len(other_z),
                    'coefs': matching_coef,
                    'odds_ratio': np.exp(matching_coef)
                })
                
    return pd.DataFrame(rows)
        
    
    

    
def assign_groups(dfDict,featuresDf,trainingDataHigh,trainingDataLow):
    '''input: groups = list of strings
              dfList = list of prs dataframes
        
        output = list of dataframes, each of which has a column for group in groups with 1 or 0 if person is assigned to group

    ''' 

    groupDfList = []
    #get exclusive main, epi and cardio
    for g,df in dfDict.items():
        try:
            df.rename(columns={'Unnamed: 0':'IID'},inplace=True)
        except KeyError:
            pass
        df.set_index('IID',inplace=True)
        #assign main exclusive
        trainingDataHigh[g] = 0
        trainingDataHigh.loc[trainingDataHigh[f'bin_{g}'] > 8,g] = 1        #assign main exclusive
        trainingDataLow[g] = 0
        trainingDataLow.loc[trainingDataLow[f'bin_{g}'] < 3,g] = 1
        df = df[featuresDf[featuresDf['model'] == g]['feature'].tolist()]
        groupDfList.append(df)
    return(groupDfList)
    
    
    

def main(phenoData,feature_scores_file,threshold):
    #########################  FILE PATHS  #################
    scoresPath = f'{phenoData}/scores'
    figPath = f'{phenoData}/figures'
    
    
    ######################## FINAL FILES  ####################
    
    featureImportanceFile = f"{scoresPath}/importantFeaturesAcrossCohortsAndTrainingData.csv"
    aucTableFile = f"{scoresPath}/performanceMetricsAcrossCohortsAndTrainingData.csv"
    
    statsTable = f"{scoresPath}/prs_pairwise_comparisons.csv"
    
    ############ IS EPI+MAIN INCLUDED IN IMPORTANT FEATURE ANALYSIS  ##############
    stats = pd.read_csv(statsTable)
    
    epiMainStats = stats[(stats['model2'] == 'scaled_prs_epi+main') & (stats['delong_p_value'] > .05)]
    otherStats = stats[(stats['model2'] != 'scaled_prs_epi+main') & (stats['delong_p_value'] > .05)]
    otherStats = otherStats[otherStats['model2'] != 'scaled_prs_all']
    
    #########################  DOWNLOAD DATAFRAMES  #################
    
    # Find PRS files
    pattern = os.path.join(f'{scoresPath}/', "*mixed*")
    tempFiles = glob.glob(pattern)
    prsFiles = [f for f in tempFiles if 'holdout' not in f and 'fromAll' not in f and 'all' not in f and 'prscr' not in f and 'covariate' not in f]
    #if epi+main not statistically different from epi or main, leave out of analysis
    if len(epiMainStats) > 0:
        prsFiles = [f for f in prsFiles if 'epi+main' not in f]
        
    if len(otherStats) > 0:
        print('other PRS models not statically significant to consider are:')
        print(otherStats)
        for row in otherStats.iterrows():
            #find the largest auc and keep that scaled_prs
            model1,auc1 = row[1]['model1'],row[1]['auc_model1']
            model2,auc2 = row[1]['model2'],row[1]['auc_model2']
            if auc2 > auc1:
                prsFiles = [f for f in prsFiles if model1.replace('scaled_prs_','') not in f]
            else:
                prsFiles = [f for f in prsFiles if model2.replace('scaled_prs_','') not in f]
                
    
    # Build dfDict dynamically
    dfDict = {}
    use_cols = []
    for file in prsFiles:
        data_type = file.split('.')[0].split('/')[-1]  # Extract just filename without path
        use_cols.append(f'bin_{data_type}')
        df = pd.read_csv(file)
        columns_to_drop = [col for col in df.columns if 'prs' in col]
        df.drop(columns=['PHENOTYPE']+columns_to_drop,inplace=True)
        dfDict[data_type] = df
        
    ################## cohort features ###################
    featuresDf = pd.read_csv(feature_scores_file)
    featuresDf = featuresDf[
        (featuresDf['coefs'] != 0) & 
        (~featuresDf['feature'].isin(['(Intercept)', 'SEX', 'age', 
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 
            'PC6', 'PC7', 'PC8', 'PC9', 'PC10']))
    ]
    
    ################## GROUPED iPRS / PRS DATAFRAME ################
    groupedDf = pd.read_csv(f'{scoresPath}/combinedPRSGroups.csv', usecols=use_cols+['IID','PHENOTYPE'])
    groupedDf.set_index('IID', inplace=True)
    
    
    
    ######################## DEFINE GROUPS DYNAMICALLY ##################
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
    
    ######################## DEFINE THE GROUPS ##################

    groups = dfDict.keys()
    dfList = dfDict.values()
    groupDfList = assign_groups(dfDict,featuresDf,trainingDataHigh,trainingDataLow)
    featuresDf.rename(columns={'model':'features_used'},inplace=True)

    trainingData = pd.concat([trainingDataHigh,trainingDataLow],ignore_index=False)
    featureImportanceFinal = pd.DataFrame()
    aucTableFinal = pd.DataFrame()
    
    for trainingDataTuple in [(trainingData,'CasesControls'),(trainingDataHigh,'HighCases'),(trainingDataLow,'LowControls')]:
        trainingData = trainingDataTuple[0]
        trainingData_str = trainingDataTuple[1]
        
        print(f'########################  TRAINING {trainingData_str}  ########################')
        print('')
        print('')
        probDf = pd.DataFrame()
        
        featureImportanceDf = pd.DataFrame()
        aucTable = pd.DataFrame()
        
        for data_type,df in dfDict.items():
            #capture auc using weighted features main, epi, epi+main, cardio features separately, so  

            feature_str = data_type        
            print(f'###################  FEATURE IMPORTANCE WITH {feature_str} FEATURES ###################')

            
            y = trainingData[groups]
            X = trainingData[['PHENOTYPE']].merge(df,left_index=True,right_index=True).drop(columns=['PHENOTYPE'])
            
            features = X.columns
            
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)
            
            # Multi-label classification using Random Forest
            model = MultiOutputClassifier(RandomForestClassifier(n_estimators=1000))
            model.fit(X_train, y_train)
            
            # Predictions
            y_pred = model.predict(X_test)
            
            # Evaluation
            report_dict = classification_report(y_test, y_pred, target_names=groups,output_dict=True)
            
            report_df = pd.DataFrame(report_dict).transpose().reset_index()

            
            for i,cl in enumerate(groups):
                probDfClass = pd.DataFrame(data=model.predict_proba(X_test)[i][:,1])
                probDfClass.columns = [f'yHat.{cl}.features.{feature_str}']

                #calculate AUC for class
                # Compute ROC curve
                fpr, tpr, _ = roc_curve(y_test[cl],model.predict_proba(X_test)[i][:,1])
                
                # Compute AUC
                roc_auc = auc(fpr, tpr)
                
                #assign auc to the report dataframe
                report_df.loc[report_df['index'] == cl,'auc'] = roc_auc
                
                print(f'AUC for {cl}: {roc_auc}')
                
                # Plot ROC curve
                plt.figure()
                plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC for {cl} using {feature_str} features (area = {roc_auc:.2f})')
                plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
                plt.xlabel('FPR')
                plt.ylabel('TPR')
                plt.title(f'AUC for class {cl} using features {feature_str}')
        #       plt.legend(loc="lower right")
                plt.savefig(f'{figPath}/{cl}.features.{feature_str}.{trainingData_str}.AUC.png')
                plt.close()
        
                print(f"Feature importance for label {cl}:")
                # Example: Analyzing feature importance for the first output
                
                explainer = fasttreeshap.TreeExplainer(model.estimators_[i], algorithm='auto',n_jobs=-1)
                shap_explainer = explainer(X_test,check_additivity=False)
                shap_values = shap_explainer.values
                
                # Visualize feature importance
        #       shap.summary_plot(shap_values[:,:,1], X_test)
                
                shap_df = pd.DataFrame(data=shap_values[:,:,1],columns=features)
                create_summary_plots(figPath,shap_values[:,:,1],X_test,i,f'{cl}.{feature_str}.{trainingData_str}.')
                
                #get top features from shap values
                topFeatures,mean_shap_values,shap_z_scores = get_top_shap_features(shap_df,data_type,threshold)
                
                plot_and_save_top_features(shap_z_scores,i,figPath,data_type)
                
                
                # Combine the Series into a DataFrame
                featureImportance = pd.concat([shap_z_scores, mean_shap_values], axis=1)
                featureImportance.columns=['z_score','mean_shap']
                featureImportance.dropna(inplace=True)
                featureImportance.reset_index(names="feature",inplace=True)
                
                #assign important feature binary to important features
                featureImportance['feature_importance'] = 0
                featureImportance.loc[featureImportance['feature'].isin(topFeatures.index.tolist()),'feature_importance'] = 1
                
                featureImportance['cohort'] = cl
                featureImportance['features_used'] = feature_str
                featureImportance['training_data_used'] = trainingData_str
                featureImportance = featureImportance.merge(featuresDf,on=['feature','features_used'],how='left')
                
                probDf = pd.concat([probDf,probDfClass],axis=1)
                
                ############# SAVE FEATURE IMPORTANCE FILE ###############
                
                if not os.path.exists(featureImportanceFile):
                    with open(featureImportanceFile,mode='w',newline='') as f:
                        featureImportance.to_csv(f,index=False)
                        f.close()
                        
                else:
                    with open(featureImportanceFile,mode='a',newline='') as f:
                        featureImportance.to_csv(f,index=False,header=False)
                        f.close()
                        
            ############# SAVE AUC REPORT TABLE ###############
                        
            
            report_df['features_used'] = feature_str
            report_df['training_data_used'] = trainingData_str
                

                        
            if not os.path.exists(aucTableFile):
                with open(aucTableFile,mode='w',newline='') as f:
                    report_df.to_csv(f,index=False)
                    f.close()
            else:
                with open(aucTableFile,mode='a',newline='') as f:
                    report_df.to_csv(f,index=False,header=False)
                    f.close()
            
        probDf.index = X_test.index
        probDf.to_csv(f"{scoresPath}/probabilitiesAcrossCohortsAndFeaturesUsedInTraining.{trainingData_str}.csv")
        
    #get the top features for the trainingDataUsed
    importantTrainindDataFeatures = pd.read_csv(f"{scoresPath}/importantFeaturesAcrossCohortsAndTrainingData.csv")
    
    filteredImportantFeaturesFile =  f"{scoresPath}/importantFeaturesAcrossCohortsAndTrainingData.Filtered.csv"   
    filteredResults = pd.DataFrame()
    
    for trainingData_str in importantTrainindDataFeatures['training_data_used'].unique():
        trainingDataFeatures = importantTrainindDataFeatures[importantTrainindDataFeatures['training_data_used'] == trainingData_str]
        feature_sets = trainingDataFeatures['features_used'].unique()
        
        # Usage
        results_df = find_opposite_sign_features(trainingDataFeatures, feature_sets,.5)
        results_df['training_data_used'] = trainingData_str
        
        # Sort by odds ratio (descending) to see strongest effects first
        results_df = results_df.sort_values('odds_ratio', ascending=False)
        
        filteredResults = pd.concat([filteredResults,results_df],ignore_index=True)
            
    filteredResults.to_csv(f"{scoresPath}/importantFeaturesAcrossCohortsAndTrainingData.Filtered.csv",index=False)


    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Calculating important features in statistically significant cohorts ....")
    parser.add_argument("--pheno_data",help="path to pheno results directory")
    parser.add_argument("--feature_scores_file",help="path to feature scores file used to calculate prs and used to find unique features for cohort")
    parser.add_argument("--threshold",type=float,default=1.99, help="Z-score threshold for identifying important features (default: 2.0)")
    
    
    args = parser.parse_args()
    
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
    feature_scores_file = args.feature_scores_file or os.environ.get("FEATURE_SCORES_FILE")
    print(f"[PYTHON] Reading features from: {feature_scores_file}")
    
    threshold = os.environ.get("THRESHOLD")
    threshold = float(threshold) if threshold else args.threshold
    
    pheno = 'celiacDisease'
    pheno_data = f'/Users/kerimulterer/prsInteractive/results/{pheno}/summedEpi'
    feature_scores_file = f'{pheno_data}/scores/featureScoresReducedFinalModel.filtered.csv'
    threshold = 1.99
    
    if not pheno_data:
        raise ValueError("You must provide a data pheno path via --pheno_data or set the PHENO_DATA environment variable.")
        
    if not feature_scores_file:
        raise ValueError("You must pass in a feature scores file with beta coefficients from association model via --feature_scores_file or set to FEATURE_SCORES_FILE environment variable")
        
    print(f"analyzing features based on shap z score : {feature_scores_file}")
    
    main(pheno_data,feature_scores_file,threshold)