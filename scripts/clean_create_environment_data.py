#!/usr/bin/env python3


def mean_center_data(df):
    '''mean center data for all features in featureList
        input : df with individuals in training set to be mean centered 
        output : df with name of column preserved but is mean centered for those features in list'''
    mean = df.mean(axis=0)
    dfMeanCentered = df - mean
    
    return(dfMeanCentered,mean)

def get_important_gene_env_features(cardioGeneticFeaturePath):
    '''input : filepath to .csv post geneticEnvironment modelling, important features post SHAP non-additive analysis
    
       output : list of env features with shap values > 2 on their own and not part of GxE interactions
                list of GxGxE interactions [[G,E],[GxG,E]...]
    '''
    
    #columns = geneticEnvFeature	genoData	Envfeature  main_E_effect
    cardioGeneticDf = pd.read_csv(cardioGeneticFeaturePath)
    mainE = cardioGeneticDf[cardioGeneticDf['main_E'] == 1]
    combinedGeneticE = cardioGeneticDf[cardioGeneticDf['main_E'] == 0]
    return(mainE,combinedGeneticE)

def combine_gene_environment(cardioFull,geneticDf,combinedGeneticE):
    pass


def scale_combine_epi_genetic_features(environFull,cardioGeneticFeaturePath,trainingDf):

    envCopy = environFull.copy()
    
    #get the genetic-env features post modelling
    cardioGeneticDf = pd.read_csv(cardioGeneticFeaturePath)
    
    
    envCopy = envCopy[envCopy.index.isin(trainingDf.index)]
    
    envCopyImputed = impute_data(envCopy)
    
    
    
    for epiCombo in epiCardioFeaturesToCombineList:
        XcenteredInteractions[','.join(epiCombo[:])] = Xcentered[epiCombo[0]] * X[epiCombo[1]]
        
    XepiScaled,scalerEpi = scale_cardio_training_data(XcenteredInteractions)
    
    #scale the epiCentered data after combining
    Xcentered,scalerEpiCentered = scale_cardio_training_data(Xcentered)
    
    Xmain = trainingCardioDf[mainCardioFeatures]
    if not Xmain.empty:
        XmainScaled,scalerMain = scale_cardio_training_data(Xmain)
        
    else:
        print("Main cardio DataFrame is empty. Scaling is skipped.")
        scalerMain = None
        
    #save data in dataframe
    #save means for epi features
    df = pd.DataFrame(data=epiMean)
    df.reset_index(inplace=True)
    df.columns = ['feature','mean']
    models = pd.DataFrame({'scalerEpi':scalerEpi,'scalerMain':scalerMain,'scalerEpiCentered':scalerEpiCentered},index=[0])
    df.to_csv(f'{dataPath}/data/cardioTrainingMeans.csv',index=False)
    models.to_csv(f'{dataPath}/data/cardioTrainingModels.csv',index=False)
        
    return(epiMean,scalerEpi,scalerMain,scalerEpiCentered)


def process_cardio_data(X,cardioDf,cardioDfFull,epiCardioFeaturesToCenter,mainCardioFeatures,epiCardioFeaturesToCombineList,dataPath,full_columns,get_dataset):
    '''transform cardioDf and geno dataset into a combined, mean centered, standardized dataset
    
    input : X dataframe(training or test, genotyped dataset with HLA region
            epiCardioFeaturesToCenter : list of cardio features involved in epistatic interactions that need to be mean centered before combining with geno data
            mainCardioFeatures : list of main cardio features not involved in epi interactions and dont need mean centering
            epiCardioFeaturesToCombineList : zipped list of cardio features and geno features in each epi pair to create new feature for model
            
    return:
            Xtransformed : dataframe combined centered and standardized dataset 
    '''	
    
    mean,scalerEpi,scalerMain,scalerEpiCentered = calculate_mean_scaler_for_cardio_training_data(cardioDfFull,epiCardioFeaturesToCenter,mainCardioFeatures,epiCardioFeaturesToCombineList,dataPath,full_columns,get_dataset)
    
    #mean center the epiCardioFeatures before combining with genotyped data for model training
    cardioMain = cardioDf[mainCardioFeatures]
    cardioEpi = cardioDf[epiCardioFeaturesToCenter]
    
    
    #mean center main	
    Xcentered = cardioEpi-mean
    
    #multiply the mean centered cardio features with geno features for seed epi pairs
    XcenteredInteractions = pd.DataFrame(index=Xcentered.index)
    
#	#get the uncentered data
#	XuncenteredInteractions = pd.DataFrame(cardioDf[epiCardioFeaturesToCenter].index)
    
    for epiCombo in epiCardioFeaturesToCombineList:
        XcenteredInteractions[','.join(epiCombo[:])] = Xcentered[epiCombo[0]] * X[epiCombo[1]]
        
#	for epiCombo in epiCardioFeaturesToCombineList:
#		XuncenteredInteractions[','.join(epiCombo[:])] = cardioDf[epiCombo[0]] * X[epiCombo[1]]
        
    # Transform the data using the fitted scaler
    scaled_data = scalerEpi.transform(XcenteredInteractions)
    
    # Create a new DataFrame with the scaled data
    XepiScaled = pd.DataFrame(scaled_data, columns=XcenteredInteractions.columns,index=XcenteredInteractions.index)
    
    # Transform epi centered data to scaled epi centered data
    scaled_data = scalerEpiCentered.transform(Xcentered)
    XcenteredScaled = pd.DataFrame(scaled_data, columns=Xcentered.columns,index=Xcentered.index)
    
    
    if not cardioMain.empty:
        if scalerMain:
            #transform Xmain
            scaled_data = scalerMain.transform(cardioMain)
            # Create a new DataFrame with the scaled data
            XmainScaled = pd.DataFrame(scaled_data, columns=cardioMain.columns,index=cardioMain.index)
        else:
            print('There was no scaler for main cardio features, in which case there were no cardio main features to scale')
            
    else:
        print("Main cardio DataFrame is empty. Scaling is skipped.")
        XmainScaled = pd.DataFrame()
        
    return(XmainScaled,XepiScaled,XcenteredScaled)

def main(resultsPath):
    
    envDf = pd.read_csv(f'{resultsPath}/participant_environment.csv')
    
    #########   GET G AND GXG FEATURES THAT ARE RANKED AND FILTERED TO COMBINE WITH E ##########
    
    trainingID = pd.read_csv(f'{resultsPath}/trainingID.txt',sep=' ')
    testID = pd.read_csv(f'{resultsPath}/testID.txt',sep=' ')
    holdoutID = pd.read_csv(f'{resultsPath}/holdoutID.txt',sep=' ')
    
    trainingEnv = envDf[envDf['IID'].isin(trainingID[0].tolist())]
    testEnv = envDf[envDf['IID'].isin(testID[0].tolist())]
    holdoutEnv = envDf[envDf['IID'].isin(holdoutID[0].tolist())]
    
    #impute missing data and mean center based on training data
    
    


    