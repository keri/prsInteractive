#!/usr/bin/env python3

import pandas as pd
import os

results_path = os.environ.get("RESULTS_PATH")
data_path = os.environ.get("DATA_PATH")



#root_dir = '/Users/kerimulterer/prsInteractive/'
hla_path = f'{results_path}/participant_hla.csv'
env_path = f'{results_path}/participant_environment.csv'
cov_path = f'{results_path}/covar.csv'
cov_path_txt = f'{results_path}/covar.txt'
withdrawals_path_txt = f'{data_path}/withdrawalsID.txt'

def convert_ID_for_csv(input_file):
	try:
		df = pd.read_csv(input_file)
		try:
			df['Participant ID'] = df['Participant ID'].str.replace('per', '', regex=False).astype(int)
		except KeyError:
			df['IID'] = df['IID'].str.replace('per', '', regex=False).astype(int)
		df.to_csv(input_file,index=False)
	except AttributeError:
		print('cannot find ',input_file)
		pass
		
def convert_ID_for_txt(input_file):
	try:
		df = pd.read_csv(input_file,sep=" ",header=False)
		try:
			df[0] = df[0].str.replace('per', '', regex=False).astype(int)
			df[1] = df[1]
		except KeyError:
			df['IID'] = df['IID'].str.replace('per', '', regex=False).astype(int)
		df.to_csv(input_file,index=False)
	except AttributeError:
		print('cannot find ',input_file)
		pass

	

for input_file in [hla_path,env_path,cov_path]:
	convert_ID_for_csv(input_file)
	
#withdrawals = pd.DataFrame([181,200,30,89])
#withdrawals.to_csv(withdrawal_path,index=False,header=None)
	
	
print('done creating IID and FID numeric values for training workflow !!')