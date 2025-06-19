#!/usr/bin/env python3

import pandas as pd
import os

root_dir = os.environ.get("PRS_INTERACTIVE_HOME")
pheno = os.environ.get("PHENO")


hla_path = f'{root_dir}/results/participant_hla.csv'
env_path = f'{root_dir}/results/participant_environment.csv'



def convert_ID_for_csv(input_file):
	try:
		df = pd.read_csv(input_file)
		df['Participant ID'] = df['Participant ID'].str.replace('per', '', regex=False).astype(int)
		df.to_csv(input_file,index=False)
	except AttributeError:
		print('cannot find ',input_file)
		pass

	

for input_file in [hla_path,env_path]:
	convert_ID_for_csv(input_file)
	
	
print('done creating IID and FID conversion file for training workflow !!')