#!/usr/bin/env python3

import pandas as pd
import numpy as np

from sklearn.metrics import roc_auc_score, make_scorer, roc_curve, auc
from sklearn.impute import SimpleImputer
import warnings
import glob
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

# Suppress the specific warning
warnings.filterwarnings("ignore", message="Bad circle positioning")


def create_feature_importance_df(model,features):
    
    feature_importances = model.feature_importances_
    
    feature_importanceDf = pd.DataFrame({'feature': features, f'feature_importance': feature_importances})
    #feature_importances = feature_importances.sort_values(by='feature_importance_gbc_mean', ascending=False)
    
    return(feature_importanceDf)

def draw_forest_plots(filePath,features,df,figPath):
    cov_features = features.tail(12)['feature'].tolist()
    
    features = features[~features['feature'].isin(cov_features)]
    features = features[features['feature'] != '(Intercept)']
    
    mainFeatures = features[features['model'] == 'main']['feature'].tolist()
    epiFeatures = features[features['model'] == 'epi']['feature'].tolist()
    epiMainFeatures = features[features['model'] == 'epi+main']['feature'].tolist()
    cardioFeatures = features[features['model'] == 'cardio']['feature'].tolist()
    
    dfPlot = df[df['model'] == 'shapley_values'][['cohort','feature','feature_importance','support']]
    
    
    #get the correct features for each cohort
    filteredPlot = dfPlot[((dfPlot['cohort'] == 'epi+main') & (dfPlot['feature'].isin(epiMainFeatures))) | ((dfPlot['cohort'] == 'main') & (dfPlot['feature'].isin(mainFeatures))) | ((dfPlot['cohort'] == 'epi') & (dfPlot['feature'].isin(epiFeatures))) | ((dfPlot['cohort'] == 'cardio') & (dfPlot['feature'].isin(cardioFeatures)))]
    
    #get the important features that have support 
    importantFeaturesDf = filteredPlot[(filteredPlot['support'] == True) & (filteredPlot['feature_importance'] > 0)]
    

    
    #get the highest feature importance for those features that are found in more than one cohort
    importantFeaturesDf.sort_values(['feature','feature_importance'],ascending=False,inplace=True)
    
    #drop duplicates, keeping the first instance
    importantFeaturesDf.drop_duplicates(subset=['feature'],keep='first',inplace=True)
    
    importantFeaturesDf.to_csv(f'{filePath}/featureImportance/groupedCohortFeatureImportanceAcrossMethods.Filtered.csv',index=False)
    importantFeaturesDf = importantFeaturesDf[importantFeaturesDf['feature_importance'] > 1]
    
    create_important_feature_forest_plot(df,importantFeaturesDf,figPath)

def calculate_important_features(X,y,allFeatures,prsPath,imageStr,figPath,balanced=1):

    #get the cross validation score for different models
    model_scores = pd.DataFrame()
    
    #get feature labels for later
    features = X.columns
    
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    Ximp = imp_mean.fit_transform(X)
    
    #the balanced dataset will have been scaled already
    if balanced == 0:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(Ximp)

    else: #this dataset is already scaled
        X_scaled = X.values
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    
    # Define AUC scoring function for multiclass
    auc_scorer = make_scorer(roc_auc_score, multi_class='ovr', needs_proba=True)
    
    # Perform cross-validation and calculate AUC for each fold
    auc_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring=auc_scorer)
    print(f'Cross-Validated AUC scores: {auc_scores}')
    print(f'Mean AUC score: {np.mean(auc_scores):.2f}')
    model_score = pd.DataFrame({'model':'random forest classifier','mean_auc':np.mean(auc_scores)},index=[0])
    model_scores = pd.concat([model_scores,model_score],ignore_index=True)
    
    # Store feature importances
    feature_importances = np.zeros(X.shape[1])
    feature_importances_relief = np.zeros(X.shape[1])
    group_probabilities_n = 0
    combined_top_features = pd.DataFrame()
    
    if balanced == 0:
        #Balance the dataset using SMOTE
        smote = SMOTE(random_state=42)
        X_resampled,y_resampled = smote.fit_resample(Ximp,y) 
    else:
        X_resampled = X.copy()
        y_resampled = y.copy()
        
        # Perform cross-validation and calculate AUC for each fold
    auc_scores = cross_val_score(rf, X_resampled, y_resampled, cv=cv, scoring=auc_scorer)
    print(f'Cross-Validated AUC scores: {auc_scores}')
    print(f'Mean AUC score: {np.mean(auc_scores):.2f}')
    model_score = pd.DataFrame({'model':'random forest classifier : balanced data','mean_auc':np.mean(auc_scores)},index=[0])
    model_scores = pd.concat([model_scores,model_score],ignore_index=True)

    
    featureImportanceDf = pd.DataFrame()
    split = 0
    for train_idx, test_idx in cv.split(X_resampled, y_resampled):
        split += 1
        X_train, X_test = X_resampled[train_idx], X_resampled[test_idx]
        y_train, y_test = y_resampled.values[train_idx], y_resampled.values[test_idx]
        
        rf.fit(X_train,y_train)

        
        ##########################  FEATURE IMPORTANCE USING SHAP VALUES  ##################
        explainer = fasttreeshap.TreeExplainer(rf, algorithm='auto',n_jobs=-1)
        shap_explainer = explainer(X_test,check_additivity=False)
        shap_values = shap_explainer.values
        #       if group_probabilities_n == 0:
        for i in range(len(rf.classes_)):
            #get the cohort for the class to calculate mean shapley values
            ind_cohort = shap_values[y_test == rf.classes_[i]]
            
            #the shape of shap_values array = (X_test.shape[0], n features, n classes)
            #get the n people x n features for each class separately
            shap_df = pd.DataFrame(data=ind_cohort.T[i].T,columns=features)
            #get top features from shap values
            topFeatures,mean_shap,z_scores_shap = get_top_features(shap_df,'topFeatures')
            featureImportance = pd.DataFrame(data=z_scores_shap,columns=['z_score'])
            featureImportance['mean_shap'] = mean_shap
            featureImportance.reset_index(inplace=True)
            featureImportance.rename(columns={'index':'feature'},inplace=True)
            featureImportance['feature_importance'] = 0
            featureImportance.loc[featureImportance['feature'].isin(topFeatures.index.tolist()),'feature_importance'] = 1
            featureImportance['cohort'] = rf.classes_[i]
            featureImportance['model'] = 'shapley_values'
            featureImportance['split'] = split
            featureImportanceDf = pd.concat([featureImportance,featureImportanceDf],ignore_index=True)

    groupedTopFeatures = featureImportanceDf.groupby(['feature','cohort','model']).mean().reset_index()

    
    ######################## FEATURE IMPORTANCE VIA RECURSIVE ELIMINATION FUNCTION ########
    #input a dataframe with feature names into the recursive function
    X_resampled_df = pd.DataFrame(data=X_resampled,columns=features)
    
    #create temp dataframe 
    featureImportance = pd.DataFrame()
    
    #run the recursive model
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    selector = RFECV(rf, min_features_to_select=100, step=5,cv=5,n_jobs=18)
    selector = selector.fit(X_resampled_df, y_resampled)
    
    #get the ranking
    ranking = selector.ranking_
    
    #get the support (True, False for each feature)
    feature_support = selector.support_
    
    #get predicted cohorts
    probs = selector.predict_proba(X_resampled)
    probDf = pd.DataFrame(data=probs)
    probDf.columns = selector.classes_
    #   probDf['cohort'] = probDf[['main','epi','epiMain']].idxmax(axis=1)
    probDf.to_csv(f'{prsPath}/recursivePredictProba.csv',index=False)
    
    # Get the class with the highest probability for each sample
    #   predicted_classes = np.argmax(probs, axis=1)
    featureImportance['feature'] = features
    featureImportance['ranking'] = ranking
    featureImportance['support'] = feature_support
    
    #merge to feature importance dataframe
    featureImportanceDf2 = featureImportanceDf.merge(featureImportance,on='feature',how="left")
    groupedTopFeatures2 = groupedTopFeatures.merge(featureImportance,on='feature',how="left")
    
    ####################### end feature importance for each cross using recursive feature elimination ########

    featureImportance.to_csv(f'{prsPath}/featureImportance/recursiveFeatureImportancePRScrGroups.csv',index=False)
    featureImportanceDf2.to_csv(f'{prsPath}/featureImportance/cohortFeatureImportanceAcrossMethodsPRScrGroups.csv',index=False)
    groupedTopFeatures2.to_csv(f'{prsPath}/featureImportance/groupedCohortFeatureImportanceAcrossMethodsPRScrGroups.csv',index=False)
#   
#   featureImportanceDf2 = pd.read_csv(f'{prsPath}/featureImportance/cohortFeatureImportanceAcrossMethodsPRScrGroups.csv')
#   groupedTopFeatures2 = pd.read_csv(f'{prsPath}/featureImportance/groupedCohortFeatureImportanceAcrossMethodsPRScrGroups.csv')
    draw_forest_plots(prsPath,allFeatures,groupedTopFeatures2,figPath)

    
def create_shap_plots(shap_values,X_train,feature_names,figPath,group):
    
    # Create the SHAP summary plot
    plt.figure(figsize=(15,10))  # Start a new figure
    shap.summary_plot(shap_values, X_train, show=False)
    
    # Save the plot as a file
    plt.savefig(f"{figPath}/shap_summary_plot.{group}.png", dpi=300, bbox_inches="tight")
    plt.close()  # Close the figure
    
def assign_groups(groupList,dfList,featuresDf,trainingDataHigh,trainingDataLow):
    '''input: groups = list of strings
              dfList = list of prs dataframes
        
        output = list of dataframes, each of which has a column for group in groups with 1 or 0 if person is assigned to group

    ''' 

    groupDfList = []
    #get exclusive main, epi and cardio
    for g,df in zip(groupList,dfList):
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
    
    
    

def main(pheno_path):
    #########################  FILE PATHS  #################
    scores_path = f'{pheno_path}/scores'
    figure_path = f'{pheno_path}/figures'

    
    ######################## FINAL FILES  ####################
    
    featureImportanceFile = f"{scores_path}/featureImportance/importantFeaturesAcrossCohortsAndTrainingData.csv"
    aucTableFile = f"{scores_path}/featureImportance/performanceMetricsAcrossCohortsAndTrainingData.csv"
    
    #########################  DOWNLOAD DATAFRAMES  #################
    
    ############## PRS DATAFAME ###############
        
    
#   mainEpiPrs = pd.read_csv(f'{prsPath}/epi+main.1539.SeparateModels.mixed.prs.csv').set_index('IID').drop(columns=['PHENOTYPE','prs','scaled_prs'])
#   mainEpiDf = pd.read_csv(f'{prsPath}/epi+main.1539.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
#   mainEpiDf.rename(columns={'Unnamed: 0':'IID'},inplace=True)
#   mainEpiDf.set_index('IID',inplace=True)
    
    if pheno == 'type2Diabetes':
        mainDf = pd.read_csv(f'{prsPath}/main.1078.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        epiDf = pd.read_csv(f'{prsPath}/epi.642.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        epiMainDf = pd.read_csv(f'{prsPath}/epi+main.1539.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        cardioDf = pd.read_csv(f'{prsPath}/cardio.221.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        groups = ['main','epi','epi+main','cardio']
        dfList = [mainDf,epiDf,epiMainDf,cardioDf]
        use_cols = ['IID','PHENOTYPE','bin_epi','bin_main','bin_cardio','bin_epi+main']
        
    else:
        mainDf = pd.read_csv(f'{prsPath}/main.372.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        epiDf = pd.read_csv(f'{prsPath}/epi.936.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        cardioDf = pd.read_csv(f'{prsPath}/cardio.207.mixed.prs.csv').drop(columns=['PHENOTYPE','prs','scaled_prs'])
        groups = ['main','epi','cardio']
        dfList = [mainDf,epiDf,cardioDf]
        use_cols = ['IID','PHENOTYPE','bin_epi','bin_main','bin_cardio']
        
    
    

#       
#   mainDf.set_index('IID',inplace=True)
#   epiDf.set_index('IID',inplace=True)
#   epiMainDf.set_index('IID',inplace=True)
##   mainEpiDf = mainDf.merge(epiDf,left_index=True,right_index=True,suffixes=['_main','_epi'])
#
#   cardioDf.set_index('IID',inplace=True)
    
    
    ################## cohort features  ###################
    
    # features within each cohort
    featuresDf = pd.read_csv(f'{filePath}/{pheno}/tanigawaSet/featureScores/featureWeights/featureScoresReducedFinalModelFilteredLD.csv')
    featuresDf = clean_feature_column_before_fixing_rscript(featuresDf)
    featuresDf = featuresDf[(featuresDf['lasso_coefs'] != 0) & (~featuresDf['feature'].isin(['(Intercept)','SEX','age','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']))]
    
    
    
    ################## GROUPED iPRS / PRS DATAFRAME WITH EPI CONTRIBUTION ANNOTATED ################
    
#   groupedDf = pd.read_csv(f'{prsPath}/combinedPRSGroups.csv',usecols=['IID','use_epi','use_main','use_epiMain','PHENOTYPE'])
    groupedDf = pd.read_csv(f'{prsPath}/combinedPRSGroupsWithcardio_main.csv',usecols=use_cols)
    groupedDf.set_index('IID',inplace=True)


    ######################## DEFINE THE GROUPS ##################
    if pheno == 'type2Diabetes':
        trainingDataHigh = groupedDf[ ((groupedDf['bin_main'] > 8) | (groupedDf['bin_epi'] > 8) | (groupedDf['bin_cardio'] > 8) | (groupedDf['bin_epi+main'] > 8)) & (groupedDf['PHENOTYPE'] == 2) ]
        trainingDataLow = groupedDf[ ((groupedDf['bin_main'] < 3) | (groupedDf['bin_epi'] < 3) | (groupedDf['bin_cardio'] < 3) | (groupedDf['bin_epi+main'] < 3)) & (groupedDf['PHENOTYPE'] == 1) ]
    else:
        trainingDataHigh = groupedDf[((groupedDf['bin_main'] > 8) | (groupedDf['bin_epi'] > 8) | (groupedDf['bin_cardio'] > 8)) & (groupedDf['PHENOTYPE'] == 2)]
        trainingDataLow = groupedDf[((groupedDf['bin_main'] < 3) | (groupedDf['bin_epi'] < 3) | (groupedDf['bin_cardio'] < 3)) & (groupedDf['PHENOTYPE'] == 1)]
    
    groupDfList = assign_groups(groups,dfList,featuresDf,trainingDataHigh,trainingDataLow)
    featuresDf.rename(columns={'model':'features_used'},inplace=True)

    trainingData = pd.concat([trainingDataHigh,trainingDataLow],ignore_index=False)
    featureImportanceFinal = pd.DataFrame()
    aucTableFinal = pd.DataFrame()
    for trainingDataTuple in [(trainingData,'CasesControls'),(trainingDataHigh,'HighCases')]:
        trainingData = trainingDataTuple[0]
        trainingData_str = trainingDataTuple[1]
        
        print(f'########################  TRAINING {trainingData_str}  ########################')
        print('')
        print('')
        probDf = pd.DataFrame()
        
        featureImportanceDf = pd.DataFrame()
        aucTable = pd.DataFrame()
        
        for j,df in enumerate(groupDfList):
            #capture auc using weighted features main, epi, epi+main, cardio features separately, so  

            feature_str = groups[j]            
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
                plt.savefig(f'{figPath}{cl}.features.{feature_str}.{trainingData_str}.AUC.png')
        
        
                print(f"Feature importance for label {cl}:")
                # Example: Analyzing feature importance for the first output
                explainer = fasttreeshap.TreeExplainer(model.estimators_[i], algorithm='auto',n_jobs=-1)
                shap_explainer = explainer(X_test,check_additivity=False)
                shap_values = shap_explainer.values
                
                # Visualize feature importance
        #       shap.summary_plot(shap_values[:,:,1], X_test)
                
                shap_df = pd.DataFrame(data=shap_values[:,:,1],columns=features)
                create_shap_plots(shap_values[:,:,1],X_test,features,figPath,f'{cl}.{feature_str}.{trainingData_str}.')
                
                
                #get top features from shap values
                topFeatures,mean_shap_values,shap_z_scores = get_top_features(shap_df,'topFeatures')
                
                create_plots(shap_df,shap_explainer,X_train,figPath,cl,f'features.{feature_str}.{trainingData_str}')
                
                plot_and_save_top_features(topFeatures,cl,f'features.{feature_str}.{trainingData_str}',figPath)
                
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
        probDf.to_csv(f"{prsPath}/featureImportance/probabilitiesAcrossCohortsAndFeaturesUsedInTraining.{trainingData_str}.csv")


            
        
        #concatenate with final Df
#       featureImportanceFinal = pd.concat([featureImportanceDf,featureImportanceFinal],ignore_index=True)
#       aucTableFinal = pd.concat([aucTable,aucTableFinal],ignore_index=True)

#   featureImportanceFinal.rename(columns={'index':'feature'},inplace=True)
#   featureImportanceFinal = featureImportanceFinal.merge(features[['feature','lasso_coefs']],on=['feature'],how='left')
#   featureImportanceFinal.to_csv(f"{prsPath}/featureImportance/importantFeaturesAcrossCohortsAndTrainingData.csv",index=False)
#   
#   aucTableFinal.to_csv(f"{prsPath}/featureImportance/performanceMetricsAcrossCohortsAndTrainingData.csv",index=False)
    
if __name__ == '__main__':
    
    pheno = sys.argv[1]
#   pheno = 'type2Diabetes'
    
    main(pheno)