#!/usr/bin/env python3

import pandas as pd
import os

root_dir = os.environ.get("PROJECT_ROOT")

training_path = f'{root_dir}/testResults/type2Diabetes/trainingCombined.raw'
test_path = f'{root_dir}/testResults/type2Diabetes/testCombined.raw'
holdout_path = f'{root_dir}/testResults/type2Diabetes/holdoutCombined.raw'
hla_path = f'{root_dir}/testResults/participant_hla.csv'
env_path = f'{root_dir}/testResults/participant_environment.csv'


# include a try/except in case this step has alr
try:
	trainingDf = pd.read_csv(training_path,sep='\s+')
	trainingDf['IID'] = trainingDf['IID'].apply( lambda x : int(x.replace('per','')))
	trainingDf.to_csv(training_path, sep=" ", index=False)
except AttributeError:
	pass

try:
	testDf = pd.read_csv(test_path,sep='\s+')
	testDf['IID'] = testDf['IID'].apply( lambda x : int(x.replace('per','')))
	testDf.to_csv(test_path, sep=" ", index=False)
except AttributeError:
	pass


try:
	holdoutDf = pd.read_csv(holdout_path,sep='\s+')
	holdoutDf['IID'] = holdoutDf['IID'].apply( lambda x : int(x.replace('per','')))
	holdoutDf.to_csv(holdout_path, sep=" ", index=False)
except AttributeError:
	pass

try:
	envDf = pd.read_csv(env_path)
	envDf['Participant ID'] = envDf['Participant ID'].apply( lambda x : int(x.replace('per','')))
	envDf.to_csv(env_path,index=False)
except AttributeError:
	pass

try:
	hlaDf = pd.read_csv(hla_path)
	hlaDf['Participant ID'] = hlaDf['Participant ID'].apply( lambda x : int(x.replace('per','')))
	hlaDf.to_csv(hla_path,index=False)
except AttributeError:
	pass

print('done !!')