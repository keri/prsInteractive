#!/usr/bin/env python3
#Keri Multerer April 2025

#simulate UK Biobank participant for phenotype definition and covariate creation

import pandas as pd
import numpy as np

participant_columns = ["Participant ID",
    "Diagnoses - main ICD10",
    "Non-cancer illness code, self-reported | Instance 0",
    "Non-cancer illness code, self-reported | Instance 1",
    "Non-cancer illness code, self-reported | Instance 2",
    "Non-cancer illness code, self-reported | Instance 3","Year of birth",
    "Sex",
    "Genetic principal components | Array 1",
    "Genetic principal components | Array 2",
    "Genetic principal components | Array 3",
    "Genetic principal components | Array 4",
    "Genetic principal components | Array 5",
    "Genetic principal components | Array 6",
    "Genetic principal components | Array 7",
    "Genetic principal components | Array 8",
    "Genetic principal components | Array 9",
    "Genetic principal components | Array 10",
    "Recommended genomic analysis exclusions",
    "Date E11 first reported (non-insulin-dependent diabetes mellitus)",
    "Source of report of E11 (non-insulin-dependent diabetes mellitus)",
    "Age at recruitment"]

covar_columns = ["IID", "FID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "age"]

##################### CREATE AND CLEAN PARTICIPANT DATA ###############

#create participant data using IID per0 ... per19999
#use simulated data .fam file
participants = pd.read_csv("~/prsInteractive/testData/sim_variant_calls/ukb22418_c22_b0_v2.fam",header=None, sep=' ')

#change 50% of column 4 to 0 (female) and 1 (male)
participants.columns = ['IID','FID','paternal','maternal','SEX','PHENOTYPE']
#change SEX column from 2 to 0
participants['SEX'] = 0

n = len(participants)
indices_to_change = np.random.choice(participants.index, size=n // 2, replace=False)

# Set those rows' values to 1
participants.loc[indices_to_change, 'SEX'] = 1

########################### PHENOTYPE DATA ##########################

phenoDf = participants[['IID','FID','PHENOTYPE']]

phenoDf.to_csv("~/prsInteractive/testData/pheno.txt",sep=' ',index=False, header=None)

###########################  CREATE COVARIATE DATA ###############

covar = participants[['IID','FID','SEX']]
#create PC columns randomly choosing values between -10 - 10
for i in range(1,11):
    covar[f'PC{i}'] = np.random.uniform(-10, 10, len(participants)).tolist()
    
covar['age'] = np.random.randint(45, 65, len(participants))

covar.to_csv("~/prsInteractive/testData/covar.txt",sep=' ',index=False)

##################  create clinical marker data  ##############

clinical_measures = {'HbA1c':[1,27,.001],'Glucose':[15,172,.01],'BMI':[15,75,.005],'M1':[1,36,.004],'M2':[1,20,.05],'M3':[90,250,.03],'M4':[10,50,.1]}

environmentalData = participants[['IID']]

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

environmentalData.to_csv("~/prsInteractive/testData/environmental.raw.csv",index=False)

################################ CREATE HLA DATA ######################

hlaData = environmentalData[['Participant ID']]
    
for i in range(1,5):
    hlaData[f'hla{i}'] = np.random.uniform(0, 2, len(hlaData)).tolist()
    
hlaData.to_csv("~/prsInteractive/testData/hlaData.raw.csv",index=False)
    