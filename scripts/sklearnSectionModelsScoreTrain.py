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
import sys
import warnings
import os
set_config(transform_output="pandas")

warnings.simplefilter(action='ignore')

# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
sys.path.append("/Users/kerimulterer/ukbiobank/")
from batchScripts.helper.download import *


def train_models(X,y,trainingPath,pheno,data_type,i):
    '''input : X from training dataset and y
    output : pickled model saved to output file'''

    print('training models .....')
    print('shape of X in training dataset = ',X.shape)
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean').fit(X)
    Ximp = imp_mean.transform(X)
    #model path
    pickle.dump(imp_mean, open(f'{trainingPath}/models/imp_mean_{data_type}_{i}.pkl', 'wb'))
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
    pickle.dump(clfNVB, open(f'{trainingPath}/models/sklearnNaiveBayes_{data_type}_{i}.pkl', 'wb'))
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
    pickle.dump(grid_search_hgb, open(f'{trainingPath}/models/sklearnGradBoostHistClassifier_{data_type}_{i}.pkl', 'wb'))
    en = time.time()
    timeHGB = (en-st)/60
    print('time if took to train HGB = ',timeHGB,' minutes')
    return(imp_mean,clfNVB,grid_search_hgb)

def score_models(X,y,featurePath,pheno,data_type,modelFile,i,imp_mean,clfNVB,clfHGB):
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
    dfSnps['log_prob_no_diabetes'] = clfNVB.feature_log_prob_[0, :]
    dfSnps['log_prob_yes_diabetes'] = clfNVB.feature_log_prob_[1, :]
    dfSnps2 = convert_log_prob_to_odds(dfSnps)


    en = time.time()
    timenvb = (en-st)/60
    print('time if took to score models = ',timenvb,' minutes')

    return(dfSnps2)



def main(pheno,data_type,start,end):

    print('starting analysis .....')

    n=3000

    ####################################################################
    #                      CREATE VARIABLES                            #
    ####################################################################
#   machinePath = '/nfs/scratch/multerke'
    # machinePath = '/nesi/nobackup/vuw03336'
    machinePath = '/Users/kerimulterer'
    trainingPath = f'{machinePath}/ukbiobank/{pheno}/tanigawaSet'
    dataPath = f'{trainingPath}/data/'
    featurePath = f'{trainingPath}/featureScores/'
    #file = 'test_main.raw'
    file = 'trainingCombined_final.raw'
    main_pathway = dataPath + file
    test_file = 'testCombined_final.raw'
    #test_file = 'test_main_holdout.raw'
    test_pathway = dataPath + test_file
    epi_pathway = dataPath + 'epiFiles/'

    # create empty dataframe to capture the scores and snps in each iteration
    models = pd.DataFrame(columns=['model','test score','balanced score','auc','matthews_corrcoef','log_loss','jaccard_score','hamming_loss','f1_score','data_type','iteration'])

#   modelFile = f'{trainingPath}/models/modelScores/sklearnModelScoresSections.csv'
#
#   if not os.path.exists(f'{modelFile}'):
#       with open(modelFile,mode='w',newline='') as f:
#           models.to_csv(f,index=False)
#           f.close()


    print(main_pathway)

    full_columns = get_columns(trainingPath)

    start = int(start)
    end = int(end)
    print(end)
    #############################################################################
    #                 GET SNPS/SNP PAIRS FOR MODEL                              #
    ############################################################################# 



    for i in range(start,end+1):
        istart = i*n
        if data_type == 'main':
            if i == 0:
                sectionSnps = full_columns[6:istart+n]
            else:
                sectionSnps = full_columns[istart:istart+n]


        else:
            epiColumns = get_epi_columns(epi_pathway)
            sectionPairs = epiColumns[istart:istart+n]
            sectionSnps = get_epi_snps(sectionPairs)


        featureScoresFile = f'{featurePath}featureScores.csv'
        #save empty dataframe on first run
        allModelFeatures = pd.DataFrame(columns=['log_prob_no_diabetes','log_prob_yes_diabetes','odds_no_diabetes','odds_yes_diabetes','model#','feature', 'model'])
        if not os.path.exists(featureScoresFile):
            with open(featureScoresFile,mode='w',newline='') as f:
                allModelFeatures.to_csv(f,index=False)
                f.close()

        print('number of snps in this section = ',len(sectionSnps))



        #########################################################################
        #                                TRAINING  / TESTING                    #
        ######################################################################### 

        #train models and pickle trained models to be used in scoring

        mainArray = get_dataset(main_pathway,sectionSnps,full_columns) #main effect snps so both pathways are the same
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
            imp_mean,clfNVB,clfHGB = train_models(Xmain,y,trainingPath,pheno,data_type,i)


            testArray = get_dataset(test_pathway,sectionSnps,full_columns)

            yTest = testArray["PHENOTYPE"]
            Xtest = testArray.drop(columns=['PHENOTYPE'])

            if data_type != 'main':
                Xtest = create_epi_df(Xtest,sectionPairs)


            print('Test array shape = ',Xtest.shape)

            allModelFeatures = score_models(Xtest,yTest,featurePath,pheno,data_type,modelFile,i,imp_mean,clfNVB,clfHGB)
            allModelFeatures['model#'] = i
            allModelFeatures['feature'] = sectionSnps
            allModelFeatures['model'] = data_type
            with open(featureScoresFile,mode='a',newline='') as f:
                allModelFeatures.to_csv(f,index=False, header=False)
                f.close()

        i += 1


if __name__ == '__main__':

    pheno = sys.argv[1]
    start = sys.argv[2]
    end = sys.argv[3]
    data_type = sys.argv[4]

    main(pheno,data_type,start,end)
    