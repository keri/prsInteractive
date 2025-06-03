# OVERVIEW

The prs Insteractive pipeline was developed using genotyped, imputed HLA, and environmental data for UK Biobank participants with a European background. The pipeline is comprised of workflows for the inclusion of gene (G), gene-gene (GxGxE), and gene-environment (GxE) interaction weights into polygenic risk (PRS) calculations for complex traits. The pipeline was developed and validated for type 2 diabetes (T2D) and celiac disease (CD) but can be applied to any trait with a ICD10 code and/or substring filter present in the Non-cancer illness code [instances 0-2] fields. Specific G, GxG, and GxGxE cohorts are further analysed to identify underlying features important to each cohorts to identify different molecular pathways driving risk within the cohorts.


# FOLDER DIRECTORY

```bash

.root (prsInteractive)
├── data
│   └── variant_calls
├── figures
├── hpc
├── results
│   └── phenotype
│       └── epiFiles
│           └── preSummaryFiles
├── scripts
│   ├── helper
│   └── test
├── testData
│   └── variant_calls
├── testResults
│   └── type2Diabetes
│       └── epiFiles
│           └── preSummaryFiles
└── workflows

19 directories

```


# INPUT FILES NEEDED 


## participant.csv 

### ancestry filters used in development dataset
 Used in genetic principal components = Yes
 Outliers for heterozygosity or missing rate = NULL
 Sex chromosome aneuploidy = NULL
 Genetic kinship or other participants = NOT ANY OF Ten or more third-degree relatives identified
 Genetic ethnic groups = IS NOT NULL
 Genetic principal components | Array 1 = IS BETWEEN -20 - 40
 Genetic principal components | Array 2 = IS BETWEEN -25 - 10
File downloaded from cohort created using cohort browser on DNA nexus platform.

### Fields listed are minimum columns needed: 


| Fields | 
|----------|
| Participant ID  |
| Diagnoses - main ICD10  |
| Non-cancer illness code, self-reported \| Instance 0  |
| Non-cancer illness code, self-reported \| Instance 1  |
| Non-cancer illness code, self-reported \| Instance 2  |
| Non-cancer illness code, self-reported \| Instance 3  |
| Sex  |
| Genetic principal components \| Array 1  |
| Genetic principal components \| Array 2  |
| Genetic principal components \| Array 3  |
| Genetic principal components \| Array 4  |
| Genetic principal components \| Array 5  |
| Genetic principal components \| Array 6  |
| Genetic principal components \| Array 7  |
| Genetic principal components \| Array 8  |
| Genetic principal components \| Array 9  |
| Genetic principal components \| Array 10  |
| Age at recruitment  |

    
## raw variant calls in bed format

| variant files separated into chromosomes c(#) | 
|----------|
| ukb{project#}_c{#}_b0_v2.bed &emsp;   ukb{project#}_c{#}_b0_v2.bim &emsp; ukb{project#}_c{#}_b0_v2.fam|

## hla_participant.csv : imputed hla loci with values range from 0 - 2

### column headings are located in: ukb_hla_v2.txt 


## participant_environment.csv

| mandatory Fields | 
|----------|
| Participant ID  |
| EHF and environmental data  |

Electronic health record data for clinical markers for blood chemistry, blood counts, and cardiometabolic features.

Features used in analysis are listed in Supplemental Table S13 of thesis.

Workflow can use as input any clinical marker which includes blood counts, blood chemistry, and lifestyle data available at initial screening and follow up visits for all participants.
    

## withdrawals.csv
 A list of eid's provided by UK Biobank of people who have opted out of research. File consists of one column with no heading or index.
  
  
  
# Output files and folder structure:

```bash 

results (root)
├── covar.txt
├── pheno
│   ├── combinedID.txt
│   ├── epiFiles
│   │   └── preSummaryFiles
│   │       ├── trainingEpi.epi.cc.1
│   │       ├── trainingEpi.epi.cc..
│   │       ├── trainingEpi.epi.cc.40
│   │       ├── trainingEpi.epi.cc.summary.1
│   │       ├── trainingEpi.epi.cc.summary..
│   │       ├── trainingEpi.epi.cc.summary.40
│   │       └── trainingEpi.log
│   ├── figures
│   ├── holdoutCombined.bed
│   ├── holdoutCombined.bim
│   ├── holdoutCombined.fam
│   ├── holdoutCombined.log
│   ├── holdoutCombined.nosex
│   ├── holdoutCombined.raw
│   ├── holdoutCombinedRaw.log
│   ├── holdoutID.txt
│   ├── merged_allChromosomes.bed
│   ├── merged_allChromosomes.bim
│   ├── merged_allChromosomes.fam
│   ├── merged_allChromosomes.log
│   ├── merged_allChromosomes.snplist
│   ├── models
│   │   ├── imp_mean_main_{i}.pkl
│   │   ├── sklearnGradBoostHistClassifier_main_{i}.pkl
│   │   └── sklearnNaiveBayes_main_{i}.pkl
│   ├── pheno_config.sh
│   ├── pheno.txt
│   ├── scores
│   │   ├── featureScores.csv
│   │   ├── sklearnModelScoresSections.csv
│   │   └── importantFeaturesPostShap.csv
│   ├── testCombined.raw
│   ├── testCombinedRaw.log
│   ├── testID.txt
│   ├── trainingCombined.raw
│   ├── trainingCombinedRaw.log
│   └── trainingID.txt
├── participant_environment.csv
└── participant_hla.csv

6 directories 

```


# WORKFLOW: 

## G, GxG, and GxGxE Analysis Overview using T2D data as an example

![PRS Pipeline Workflow](READMEfigures/combinedGWASWorkflow.png)


# Important Underlying Feature Analysis using T2D data as an example

![Important Feature Pipeline Workflow](READMEfigures/importantFeatureWorkflowSHAP.png)


#   STEPS IN ANALYSIS 

## 1) check to see if a conda environment has been created in /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
- if ukb_env is not present, create ukb_env in root directory prsInteractive:
  Create environment from environment.yml file
  This will create a conda "ukb_env" folder with all of the dependencies in the prsInteractive/ directory

  - If on the raapoi hpc path/to/prsInteractive, 
  
  - In an interactive session run:
  
  ```bash 
  
  $ cd /path/to/directory/prsInteractive
  $ module load Miniconda3/4.9.2
  $ source $(conda info --base)/etc/profile.d/conda.sh 
  $ conda env create --prefix ./ukb_env -f environment.yml
  
  ```
  
  
## 2) with conda env "ukb_env" present run the workflow:

  * pheno = phenotype spelled in camel font and no spaces (i.e. type2Diabetes, myocardialInfarction)
  * icd10 code = substring present in the UK Biobank data
  * pheno substring = will be exact spelling found in UKB data to check for if icd10 not present. this will have spaces so will need to wrap in " "
  * n = number of cores to pass to the epistatic analysis, with 40 cores being the norm and will take approximately 48 hours
  
  #run command lines:
  
  ```bash 
  
  $ cd nfs/scratch/projects/ukbiobank/prsInteractive/hpc
  
  $ sbatch run_data_cleaning_workflow_submit.sh {pheno} {icd10 code}  {"sub string"} {n}
  #i.e. sbatch run_data_cleaning_workflow_submit.sh myocardialInfarction I21 "myocardial infarction" 40
  
  ````
    
  #### sbatch run_data_cleaning_workflow_submit.sh has the following scripts
  
  ```bash
  
  # create phenotype data and train test split IDs
  $ python "${SCRIPTS_DIR}/create_pheno_train_test_split.py"
  
  # create hla (and environmental data files?)
  $ python "${SCRIPTS_DIR}/clean_environment_hla_covar_data.py"
  
  # Run the variant call cleaning
  $ bash "${SCRIPTS_DIR}/plink_clean_variant_calls.sh"
  
  # merge separate chromosome files into one
  $ bash "${SCRIPTS_DIR}/merge_chromosomes.sh"
  
  # run epistatic interaction analysis
  $ sbatch multiprocessing_fast_epistasis_submit.sh
  
  # start the queues of batch runs for single SNPs
  $ sbatch run_model_model_batches_submit.sh {pheno} main 
  
  ```
  
## 3) After epistatic analysis is complete, run the batch models with 
  
  pheno = name of folder created and entered in run_data_cleaning_workflow_submit.sh
  data_type = main (if running single SNPs) or epi (if running with epi-pairs created from epistatic analysis
  
  input files include:
  * $PHENO_FOLDER/trainingCombined.raw
  * $PHENO_FOLDER/testCombined.raw
  * $PHENO_FOLDER/epiFiles/trainingCombinedEpi.epi.cc.summary 
  
  output files include: 
  * $PHENO_PATH/scores/featureScores.csv
  * $PHENO_PATH/scores/importantFeaturesPostShap.csv
  * $PHENO_PATH/score/sklearnModelScoresSection.csv
  * $PHENO_PATH/pheno_config.sh
  
  
  batch script generates a number of hpc jobs running 5 models with of 3K epi features each
  
  ```bash 
  
  $ cd $ROOT_DIRECTORY/hpc/
  
  #scripts removes redundant epi pairs and creates a filtered summary file
  #each hpc job takes approximately 6 hours using 800GB of RAM and 50 cpus
  $ sbatch run_model_epi_models.sh pheno
  
  ```
  
  ### When the models are completed: 
  
  importantFeaturesPostShap.csv is used to run the final association models
  
  Before modelling: The gene-environment datasets for training, test, holdout datasets are created with steps:
  * impute and mean-center training data, use imputation model and mean to imput and mean-center test and holdout data
  * combine datasets and scaled training data set, use trained scaler from training data to scaled test and holdout data
  
  
  
  ```bash 
  
  $ cd $ROOT_DIRECTORY/hpc/
  
  #scripts removes redundant epi pairs and creates a filtered summary file
  #each hpc job takes approximately 6 hours using 800GB of RAM and 50 cpus
  $ sbatch run_model_epi_models.sh pheno
  
  ```
  
  # Running analysis with WDL workflow
  
  ### These instructions use the cromwell backend to compile .wdl workflow
  ###
  ### cromwell specific features include a cromwell.config file and command line instruction to start run
  ###
  
  ### Input
  
  All inputs are required for from the root directory and annotated in the pipelineInputs.json file
  Phenotype specific inputs are updated in first step
  ###
  
  ### Download cromwell 
  
  Instructions are found here: [cromwell download] (https://cromwell.readthedocs.io/en/latest/tutorials/FiveMinuteIntro/)
  
  
  #### Step 1: cd into root directory
  
  ```
  cd /prsInteractive
  ./envSetUp.sh <phenotype> <icd10> <phenotype string used in search> <n cores to run epistatic analysis>
  
  ```
  
  
  



  

