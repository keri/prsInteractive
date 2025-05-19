#Keri Multerer Dec 2023
#Models generated via 5K sections

import pandas as pd
import numpy as np
import time
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.naive_bayes import ComplementNB
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score, matthews_corrcoef, log_loss, jaccard_score, hamming_loss
from sklearn.impute import SimpleImputer
from sklearn import set_config
import pickle
import csv
import warnings
import os
import argparse
set_config(transform_output="pandas")

warnings.simplefilter(action='ignore')

# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
#sys.path.append("/Users/kerimulterer/ukbiobank/")
from helper.download import get_dataset, get_epi_columns, get_columns
from helper.calculate_shap_values import *
from helper.data_wrangling import *



def convert_log_prob_to_odds(df):
    '''input df columns = ['feature', 'model', 'log_prob_no_pheno', 'log_prob_yes_pheno',
       'risk_ratio_no_pheno', 'risk_ratio_yes_pheno', 'lasso_coef',
       'chr_section']
       output : columns added ['odds_yes_pheno','odds_no_pheno']'''
    #convert log probabilities to probabilities
    prob_no_pheno = np.exp(-(df['log_prob_no_pheno']))
    prob_yes_pheno = np.exp(-(df['log_prob_yes_pheno']))
    
    df['odds_yes_pheno'] = (df['log_prob_yes_pheno']/df['log_prob_no_pheno']) / (prob_yes_pheno/prob_no_pheno)
    df['odds_no_pheno'] = (df['log_prob_no_pheno']/df['log_prob_yes_pheno']) / (prob_no_pheno/prob_yes_pheno)
    
    df.sort_values(['odds_yes_pheno'],ascending=False,inplace=True)
    return(df)


def train_models(X,y,modelPath,pheno,data_type,i):
    '''input : X from training dataset and y
    output : pickled model saved to output file'''

    print('training models .....')
    print('shape of X in training dataset = ',X.shape)
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean').fit(X)
    Ximp = imp_mean.transform(X)
    #model path
    pickle.dump(imp_mean, open(f'{modelPath}/imp_mean_{data_type}_{i}.pkl', 'wb'))
    ##########################
    # Lasso regression
    ##########################

    # st = time.time()
    # #clfLasso = LogisticRegressionCV(cv=10, random_state=0,penalty='elasticnet',l1_ratios=np.arange(0, 1, 0.1),solver='saga',n_jobs=120).fit(Ximp, y)
    # clfLasso = LogisticRegressionCV(cv=10, random_state=0,penalty='l1',solver='liblinear',n_jobs=120).fit(Ximp, y)
    # en = time.time()
    # timeLasso = (en-st)/60
    # pickle.dump(clfLasso, open(f'{trainingPath}models/sklearnLasso_{data_type}_{i}V2.pkl', 'wb'))
    # print('time if took to train lasso = ',timeLasso,' minutes')

    ##########################
    # Naive Bayes
    ##########################

    st = time.time()
    clfNVB = ComplementNB().fit(Ximp,y)
    en = time.time()
    timenvb = (en-st)/60
    pickle.dump(clfNVB, open(f'{modelPath}/sklearnNaiveBayes_{data_type}_{i}.pkl', 'wb'))
    print('time if took to train naive bayes = ',timenvb,' minutes')

    ##########################
    # Gradient boosted Hist classifier
    ##########################
    st = time.time()
    parameters_hgb = [{'max_iter':[1000,1500,2000],'learning_rate':[.001,.01,.1,1],'l2_regularization': [0,.1,.5]}]
    # clfHGB = HistGradientBoostingClassifier(early_stopping='auto',cv=10,l2_regularization=0, learning_rate=0.001,
    #                            max_depth=25, max_iter=1000, scoring='f1_micro').fit(X, y)
    clfHGB = HistGradientBoostingClassifier(early_stopping='auto')
    grid_search_hgb = GridSearchCV(estimator=clfHGB,param_grid=parameters_hgb,scoring='roc_auc',cv=3,n_jobs=50)
    grid_search_hgb.fit(X.to_numpy(),y.to_numpy())
    pickle.dump(grid_search_hgb, open(f'{modelPath}/sklearnGradBoostHistClassifier_{data_type}_{i}.pkl', 'wb'))
    en = time.time()
    timeHGB = (en-st)/60
    print('time if took to train HGB = ',timeHGB,' minutes')
    return(imp_mean,clfNVB,grid_search_hgb)

def score_models(X,y,pheno,data_type,modelFile,i,imp_mean,clfNVB,clfHGB,figPath): 
    '''load pickled models and score with test set'''
    print('scoring models .....')

    st = time.time()

    #get the feature names for model

    Ximp = imp_mean.transform(X)

    #########################
    # Gradient Boosted Hist 
    #########################
    #model = pickle.load(open(f'{trainingPath}/models/sklearnGradBoostHistClassifier_{data_type}_{i}', 'rb'))
    #clfHGB = pickle.load(model)
    score = clfHGB.score(X.to_numpy(), y.to_numpy())
    yHat = clfHGB.predict(X.to_numpy())
    balanced_score = balanced_accuracy_score(y.to_numpy(),yHat)
    auc = roc_auc_score(y.to_numpy(), clfHGB.predict_proba(X.to_numpy())[:, 1])
    mcc = matthews_corrcoef(y.to_numpy(),yHat)
    logloss = log_loss(y.to_numpy(),clfHGB.predict_proba(X.to_numpy())[:, 1])
    jscore = jaccard_score(y.to_numpy(),yHat)
    hloss = hamming_loss(y.to_numpy(),yHat)
    f1score = f1_score(y.to_numpy(),yHat)
    fields=['gradient boosted classifier',score,balanced_score,auc,mcc,logloss,jscore,hloss,f1score,data_type,i]

    with open(modelFile,mode='a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        f.close()

#     except FileNotFoundError:
#         pass


    ##########################
    # Naive Bayes
    ##########################

#     try:

#    model = pickle.load(open(f'{trainingPath}/models/sklearnNaiveBayes_{data_type}_{i}.pkl', 'rb'))
#    clfNVB = pickle.load(model)
    if auc > .51:
        rank_features = True
    
    score = clfNVB.score(Ximp, y)
    yHat = clfNVB.predict(Ximp)
    balanced_score = balanced_accuracy_score(y,yHat)
    auc = roc_auc_score(y, clfNVB.predict_proba(Ximp)[:, 1])
    mcc = matthews_corrcoef(y,yHat)
    logloss = log_loss(y,clfNVB.predict_proba(Ximp)[:, 1])
    jscore = jaccard_score(y,yHat)
    hloss = hamming_loss(y,yHat)
    f1score = f1_score(y,yHat)
    fields=['naive bayes',score,balanced_score,auc,mcc,logloss,jscore,hloss,f1score,data_type,i]

    with open(modelFile,mode='a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        f.close()


    dfSnps = pd.DataFrame()
    dfSnps['log_prob_no_pheno'] = clfNVB.feature_log_prob_[0, :]
    dfSnps['log_prob_yes_pheno'] = clfNVB.feature_log_prob_[1, :]
    dfSnps2 = convert_log_prob_to_odds(dfSnps)


    en = time.time()
    timenvb = (en-st)/60
    print('time if took to score models = ',timenvb,' minutes')
    
    if rank_features or auc > .51:
        topFeatures,featuresZscores = calculate_plot_shap_values(clfHGB,X,y,i,figPath,data_type)
    
    
    return(dfSnps2,topFeatures)



def main(pheno,pheno_path,training_path,test_path,epi_path,data_type,start,end):

    print('starting analysis .....')

    n=3000

    ####################################################################
    #                      CREATE VARIABLES                            #
    ####################################################################
    
    scoresPath = f'{pheno_path}/scores'
    modelPath = f'{pheno_path}/models'
    figPath = f'{pheno_path}/figures'

    # create empty dataframe to capture the scores and snps in each iteration
    models = pd.DataFrame(columns=['model','test score','balanced score','auc','matthews_corrcoef','log_loss','jaccard_score','hamming_loss','f1_score','data_type','iteration'])

    modelFile = f'{scoresPath}/sklearnModelScoresSections.csv'

    if not os.path.exists(f'{modelFile}'):
        print(f'creating model scores file : {modelFile} ...')
        with open(modelFile,mode='w',newline='') as f:
            models.to_csv(f,index=False)
            f.close()


    print(pheno_path)

    full_columns = get_columns(pheno_path)

    start = int(start)
    end = int(end)
    print(end)
    #############################################################################
    #                 GET SNPS/SNP PAIRS FOR MODEL                              #
    ############################################################################# 



    for i in range(start,end+1):
        istart = (i-1)*n
        if data_type == 'main':
#           if i == 0:
#               sectionSnps = full_columns[:istart+n]
#           else:
            sectionSnps = full_columns[istart:istart+n]


        else:
            print('epi_path to get filtered epi pairs ..',epi_path)
            epiColumns = get_epi_columns(epi_path)
            sectionPairs = epiColumns[istart:istart+n]
            sectionSnps = get_epi_snps(sectionPairs)


        featureScoresFile = f'{scoresPath}/featureScores.csv'
        importantFeaturesFile = f'{scoresPath}/importantFeaturesPostShap.csv'
        
        

        if not os.path.exists(featureScoresFile):
            #save empty dataframe on first run
            allModelFeatures = pd.DataFrame(columns=['log_prob_no_pheno','log_prob_yes_pheno','odds_no_pheno','odds_yes_pheno','model#','feature', 'model'])
            
            print(f'creating features scores file : {featureScoresFile} ...')
            
            with open(featureScoresFile,mode='w',newline='') as f:
                allModelFeatures.to_csv(f,index=False)
                f.close()

        print('number of snps in this section = ',len(sectionSnps))
        
        ### save empty important features df
        if not os.path.exists(importantFeaturesFile):
            #save empty dataframe on first run
            importantFeaturesShap = pd.DataFrame(columns=['feature','shap_zscore','data_type'])
            
            print(f'creating important features file : {importantFeaturesFile} ...')
            
            with open(importantFeaturesFile,mode='w',newline='') as f:
                importantFeaturesShap.to_csv(f,index=False)
                f.close()
                
        
        #########################################################################
        #                                TRAINING  / TESTING                    #
        ######################################################################### 

        #train models and pickle trained models to be used in scoring

        mainArray = get_dataset(training_path,sectionSnps,full_columns) #main effect snps so both pathways are the same
    #     the first columns will be IID, PHENOTYPE
        y = mainArray['PHENOTYPE']
        Xmain = mainArray.drop(columns=["PHENOTYPE"])

        if data_type != 'main':
            Xmain = create_epi_df(Xmain,sectionPairs)


        print('Xmain array = ',Xmain.shape)

        print('training section model ....')
    #     train all snps
    #     full columns have [IID,FID,FAM,DISTANCE,SEX,PHENOTYPE]


        if Xmain.shape[1] > 0:
            imp_mean,clfNVB,clfHGB = train_models(Xmain,y,modelPath,pheno,data_type,i)


            testArray = get_dataset(test_path,sectionSnps,full_columns)

            yTest = testArray["PHENOTYPE"]
            Xtest = testArray.drop(columns=['PHENOTYPE'])

            if data_type != 'main':
                Xtest = create_epi_df(Xtest,sectionPairs)


            print('Test array shape = ',Xtest.shape)

            allModelFeatures,topFeatures = score_models(Xtest,yTest,pheno,data_type,modelFile,i,imp_mean,clfNVB,clfHGB,figPath,)
            allModelFeatures['model#'] = i
            allModelFeatures['feature'] = sectionSnps
            allModelFeatures['model'] = data_type
            
            if topFeatures.empty:
                pass
            else:
                topFeatures = topFeatures.reset_index()
                topFeatures['data_type'] = data_type

                with open(importantFeaturesFile,mode='a',newline='') as f:
                    topFeatures.to_csv(f,index=False, header=False)
                    f.close()
                    
            with open(featureScoresFile,mode='a',newline='') as f:
                allModelFeatures.to_csv(f,index=False, header=False)
                f.close()

        i += 1


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="running models for wrapper...")
    parser.add_argument("--pheno_folder", help="Path to the input pheno folder")
    parser.add_argument("--training_file", help="data file of training data")
    parser.add_argument("--epi_file", help="epi file with epi pairs to analyze")
    parser.add_argument("--test_file", help="data file of test data")
    parser.add_argument("--data_type", help="data type to analyze")
    parser.add_argument("--start", help="start job")
    parser.add_argument("--end", help="end job")
    parser.add_argument("--pheno", help="Phenotype to analyze")

    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
#   pheno_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes'
    pheno_path = args.pheno_folder or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")
    
#   pheno = 'type2Diabetes'
    pheno = args.pheno or os.environ.get("PHENO")
    print(f"[PYTHON] Phenotype : {pheno}")
    
#   data_type = 'epi'
    data_type = args.data_type or os.environ.get("DATA_TYPE")
    print(f"data type : {data_type}")
    
#   training_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes/trainingCombined.raw'
    training_path = args.training_file or os.environ.get("TRAINING_PATH")
    print(f"training file : {training_path}")
    
#   test_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes/testCombined.raw'
    test_path = args.test_file or os.environ.get("TEST_PATH")
    print(f"test file : {test_path}")
    
    
    
    if data_type == 'epi':
        epi_path = args.epi_file or os.environ.get("EPI_PATH")
#       epi_path = '/Users/kerimulterer/prsInteractive/testResults/type2Diabetes/epiFiles/trainingCombinedEpi.epi.cc.summary'
        if not epi_path:
            raise ValueError("You must provide a data type code via --epi_path or set the EPI_PATH environment variable.")
        print(f"epi path : {epi_path}")
    else:
        epi_path = 'None'
    
#   start = 1
    start = args.start or os.environ.get("START")
    print(f"start : {start}")
    
#   end = 1
    end = args.end or os.environ.get("END")
    print(f"end : {end}")
    
    
    
    if not pheno_path:
        raise ValueError("You must provide a data pheno path via --pheno_folder or set the PHENO_PATH environment variable.")
        
    if not pheno:
        raise ValueError("You must provide a phenotype via --pheno or set the PHENO environment variable.")
        
    if not data_type:
        raise ValueError("You must provide a data type code via --data_type or set the DATA_TYPE environment variable.")
        
    if not training_path:
        raise ValueError("You must provide a data type code via --training_path or set the TRAINING_PATH environment variable.")
        
    if not test_path:
        raise ValueError("You must provide a data type code via --test_path or set the TEST_PATH environment variable.")
        
    if not start:
        raise ValueError("You must provide a data type code via --start or set the START environment variable.")
        
    if not end:
        raise ValueError("You must provide a data type code via --end or set the END environment variable.")
        
    
    #check to see that scores, models, and figures folder has been created
    dir_path = f"{pheno_path}/models"
    os.makedirs(dir_path, exist_ok=True)
    
    dir_path = f"{pheno_path}/scores"
    os.makedirs(dir_path, exist_ok=True)
    
    dir_path = f"{pheno_path}/figures"
    os.makedirs(dir_path, exist_ok=True)

        
    main(pheno,pheno_path,training_path,test_path,epi_path,data_type,start,end)
    