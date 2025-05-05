#!/usr/bin/env python3
#Keri Multerer April 2025

#simulate UK Biobank participant for phenotype definition and covariate creation

import pandas as pd
import numpy as np
import os

# Prefer command-line input if provided; fallback to env var
data_path = os.environ.get("DATA_PATH")
print(f"[PYTHON] Reading from: {data_path}")

# Prefer command-line input if provided; fallback to env var
results_path = os.environ.get("RESULTS_PATH")
print(f"[PYTHON] writing to: {results_path}")

pheno_path = os.environ.get("PHENO_PATH")
print(f"[PYTHON] Reading from: {pheno_path}")

pheno = os.environ.get("PHENO")
print(f"[PYTHON] Phenotype : {pheno}")

#participant_columns = ["Participant ID",
#    "Diagnoses - main ICD10",
#    "Non-cancer illness code, self-reported | Instance 0",
#    "Non-cancer illness code, self-reported | Instance 1",
#    "Non-cancer illness code, self-reported | Instance 2",
#    "Non-cancer illness code, self-reported | Instance 3",
#    "Sex",
#    "Genetic principal components | Array 1",
#    "Genetic principal components | Array 2",
#    "Genetic principal components | Array 3",
#    "Genetic principal components | Array 4",
#    "Genetic principal components | Array 5",
#    "Genetic principal components | Array 6",
#    "Genetic principal components | Array 7",
#    "Genetic principal components | Array 8",
#    "Genetic principal components | Array 9",
#    "Genetic principal components | Array 10",
#    "Date E11 first reported (non-insulin-dependent diabetes mellitus)",
#    "Source of report of E11 (non-insulin-dependent diabetes mellitus)",
#    "Age at recruitment"]

covar_columns = ["IID", "FID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "age"]

##################### CREATE AND CLEAN COVARIATE DATA ###############

#create participant data using IID per0 ... per19999
#use simulated data .fam file
covar = pd.read_csv(f"{data_path}/variant_calls/ukb22418_c22_b0_v2.fam",header=None, sep=' ')

#change 50% of column 4 to 0 (female) and 1 (male)
covar.columns = ['IID','FID','paternal','maternal','SEX','PHENOTYPE']
#change SEX column from 2 to 0
covar['SEX'] = 0

n = len(covar)
indices_to_change = np.random.choice(covar.index, size=n // 2, replace=False)

# Set those rows' values to 1
covar.loc[indices_to_change, 'SEX'] = 1

########################### PHENOTYPE DATA ##########################

phenoDf = covar[['IID','FID','PHENOTYPE']]

phenoDf.to_csv(f"{results_path}/pheno.txt",sep=' ',index=False, header=None)



###########################  CREATE COVARIATE DATA ###############

covar = covar[['IID','FID','SEX']]
#create PC columns randomly choosing values between -10 - 10
for i in range(1,11):
    covar[f'PC{i}'] = np.random.uniform(-10, 10, len(covar)).tolist()
    
covar['age'] = np.random.randint(45, 65, len(covar))

covar.to_csv(f"{results_path}/covar.txt",sep=' ',index=False)

##################  create clinical marker data  ##############

clinical_measures = {'HbA1c':[1,27,.001],'Glucose':[15,172,.01],'BMI':[15,75,.005],'M1':[1,36,.004],'M2':[1,20,.05],'M3':[90,250,.03],'M4':[10,50,.1]}

environmentalData = covar[['IID']]

#change IID to Participant ID which is what the scripts were written for
environmentalData.columns = ['Participant ID']

for measure,thresholds in clinical_measures.items():
    environmentalData[measure] = np.random.randint(thresholds[0], thresholds[1], len(environmentalData))
    
    # Percentage of rows to set as NaN (e.g., 30%)
    percentage = thresholds[2]
    n_nan = int(len(environmentalData) * percentage)
    
    # Randomly choose which rows to set to NaN
    nan_indices = np.random.choice(environmentalData.index, size=n_nan, replace=False)
    
    # Set those rows in 'col' to NaN
    environmentalData.loc[nan_indices, measure] = np.nan

environmentalData.to_csv(f"{data_path}/participant_environment.csv",index=False)

################################ CREATE HLA DATA ######################

hlaData = environmentalData[['Participant ID']]
    
for i in range(1,5):
    hlaData[f'hla{i}'] = np.random.uniform(0, 2, len(hlaData)).tolist()
    
hlaData.to_csv(f"{data_path}/participant_hla.csv",index=False)

###################### CREATE PARTICIPANT DATA ########################

participants = covar.drop(columns=['FID','SEX','age']).rename(columns={'IID':'Participant ID'})

participants.columns = [col.replace('PC','Genetic principal components | Array ') for col in participants.columns]

# set the columns 
#"Diagnoses - main ICD10",
#"Non-cancer illness code, self-reported | Instance 0",
#"Non-cancer illness code, self-reported | Instance 1",
#"Non-cancer illness code, self-reported | Instance 2",
#"Non-cancer illness code, self-reported | Instance 3","Year of birth",
#"Sex",
#"age"


participants["Age at recruitment"] = covar['age']
participants["Sex"] = covar['SEX'].apply(lambda x : 'Female' if x == 1 else 'Male')

# Get all indices where phenotype == 2
phenotype_2_indices = phenoDf[phenoDf['PHENOTYPE'] == 2].index


# Create new column
participants['Non-cancer illness code, self-reported | Instance 0'] = ''
participants['Non-cancer illness code, self-reported | Instance 1'] = ''
participants['Non-cancer illness code, self-reported | Instance 2'] = ''
participants['Non-cancer illness code, self-reported | Instance 3'] = ''
participants['Diagnoses - main ICD10'] = ''
participants['Date E11 first reported (non-insulin-dependent diabetes mellitus)'] = ''
participants['Source of report of E11 (non-insulin-dependent diabetes mellitus)'] = ''


# Assign 'type 2 diabetes' to the selected to % of intances
selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.1), replace=False)
participants.loc[selected_indices, 'Non-cancer illness code, self-reported | Instance 0'] = 'type 2 diabetes'

selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.05), replace=False)
participants.loc[selected_indices, 'Non-cancer illness code, self-reported | Instance 1'] = 'type 2 diabetes'

selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.02), replace=False)
participants.loc[selected_indices, 'Non-cancer illness code, self-reported | Instance 2'] = 'type 2 diabetes'

selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.01), replace=False)
participants.loc[selected_indices, 'Non-cancer illness code, self-reported | Instance 3'] = 'type 2 diabetes'

selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.7), replace=False)
participants.loc[selected_indices, 'Diagnoses - main ICD10'] = 'E11'

selected_indices = np.random.choice(phenotype_2_indices, size=int(len(phenotype_2_indices) * 0.3), replace=False)
participants.loc[selected_indices, 'Source of report of E11 (non-insulin-dependent diabetes mellitus)'] = 'not empty'

#assign the remaining cases that aren't assigned a value as 'not empty'
columns_to_check = [ "Diagnoses - main ICD10",
    "Non-cancer illness code, self-reported | Instance 0",
    "Non-cancer illness code, self-reported | Instance 1",
    "Non-cancer illness code, self-reported | Instance 2",
    "Non-cancer illness code, self-reported | Instance 3",
    "Source of report of E11 (non-insulin-dependent diabetes mellitus)",
]

# Randomly select 1% of phenotype_indices
n_rows = int(0.005 * len(phenotype_2_indices))
selected_rows = np.random.choice(phenotype_2_indices, size=n_rows, replace=False)

# Set all values in the selected rows and columns to ''
participants.loc[selected_rows, columns_to_check] = ''


# Assign value based on the mask
participants.loc[selected_rows, 'Date E11 first reported (non-insulin-dependent diabetes mellitus)'] = 'not empty'

participants.to_csv(f"{data_path}/participant.csv",index=False)