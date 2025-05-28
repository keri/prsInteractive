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
	trainingDf['IID'] = trainingDf['IID'].str.replace('per', '', regex=False).astype(int)
	trainingDf['FID'] = trainingDf['FID'].str.replace('per', '', regex=False).astype(int)
	trainingDf.to_csv(training_path, sep=" ", index=False)
except AttributeError:
	pass

try:
	testDf = pd.read_csv(test_path,sep='\s+')
	testDf['IID'] = testDf['IID'].str.replace('per', '', regex=False).astype(int)
	testDf['FID'] = testDf['FID'].str.replace('per', '', regex=False).astype(int)
	testDf.to_csv(test_path, sep=" ", index=False)
except AttributeError:
	pass


try:
	holdoutDf = pd.read_csv(holdout_path,sep='\s+')
	holdoutDf['IID'] = holdoutDf['IID'].str.replace('per', '', regex=False).astype(int)
	holdoutDf['FID'] = holdoutDf['FID'].str.replace('per', '', regex=False).astype(int)
	holdoutDf.to_csv(holdout_path, sep=" ", index=False)
except AttributeError:
	pass

try:
	envDf = pd.read_csv(env_path)
	envDf['Participant ID'] = envDf['Participant ID'].str.replace('per', '', regex=False).astype(int)
	envDf.to_csv(env_path,index=False)
except AttributeError:
	pass

try:
	hlaDf = pd.read_csv(hla_path)
	hlaDf['Participant ID'] = hlaDf['Participant ID'].str.replace('per', '', regex=False).astype(int)
	hlaDf.to_csv(hla_path,index=False)
except AttributeError:
	pass

print('done converting IID and FID to for training workflow !!')