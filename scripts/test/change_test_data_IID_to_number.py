#!/usr/bin/env python3

import pandas as pd
import os

root_dir = os.environ.get("PRS_INTERACTIVE_HOME")


#root_dir = '/Users/kerimulterer/prsInteractive/'
hla_path = f'{root_dir}/results/participant_hla.csv'
env_path = f'{root_dir}/results/participant_environment.csv'
cov_path = f'{root_dir}/results/covar.csv'


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

	

for input_file in [hla_path,env_path,cov_path]:
	convert_ID_for_csv(input_file)
	
	
print('done creating IID and FID numeric values for training workflow !!')