# OVERVIEW

The prs Insteractive pipeline was developed using genotyped, imputed HLA, and environmental data for UK Biobank participants with a European background. The pipeline is comprised of workflows for the inclusion of gene (G), gene-gene (GxGxE), and gene-environment (GxE) interaction weights into polygenic risk (PRS) calculations for complex traits. The pipeline was developed and validated for type 2 diabetes (T2D) and celiac disease (CD) but can be applied to any trait with a ICD10 code and/or substring filter present in the Non-cancer illness code [instances 0-2] fields. Specific G, GxG, and GxGxE cohorts are further analysed to identify underlying features important to each cohorts to identify different molecular pathways driving risk within the cohorts.


# FOLDER DIRECTORY

```bash

.
├── data
│   └── variant_calls
├── hpc
├── results
│   ├── epiFiles
│   │   └── preSummaryFiles
│   ├── figures
│   └── models
├── scripts
│   ├── helper
│   └── test
├── testData
│   └── variant_calls
└── workflows

15 directories

```


# INPUT FILES IN PRS/DATA/


hla_participant.csv #imputed HLA data
participant.csv #created on DNA nexus for study cohort
participant_environment.csv #created on the DNA nexus with environment variables
ukb_hla_v2.txt #downloaded from UKBiobank showcase
variant_calls/#bed/bim/bam files for chromosomes
withdrawals.csv (used to remove withdrawals in download)
withdrawalsID.txt (used to remove withdrawals in plink steps)

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

    
## variants/raw variant calls in bed format

| variant files separated into chromosomes c(#) | 
|----------|
| ukb{project#}_c{#}_b0_v2.bed &emsp;   ukb{project#}_c{#}_b0_v2.bim &emsp; ukb{project#}_c{#}_b0_v2.fam|

## hla_participant.csv : imputed hla loci with values range from 0 - 2

### column headings are located in: ukb_hla_v2.txt 


## participant_environment.csv

### Workflow can use as input any clinical marker which includes blood counts, blood chemistry, and lifestyle data available at initial screening and follow up visits for all participants.

| mandatory Fields | 
|----------|
| Participant ID  |
| EHF and environmental data  |

- Electronic health record data for clinical markers for blood chemistry, blood counts, and cardiometabolic features.

- Features used in analysis are listed in Supplemental Table S13 of thesis.


## withdrawals.csv
### A list of eid's provided by UK Biobank of people who have opted out of research. File consists of one column with no heading or index.
  

# OUTPUT FILES AND DIRECTORIES
  
# Output files and folder structure:

```
.
.pheno/
├── finalModel.ld
├── finalModel.log
├── finalModel.prune.in
├── finalModel.prune.out
├── finalModel.tags.list
├── finalModelLDSnps.txt
├── geneEnvironmentHoldout.csv
├── geneEnvironmentTest.csv
├── geneEnvironmentTraining.csv
├── holdoutCombined.bed
├── holdoutCombined.bim
├── holdoutCombined.fam
├── holdoutCombined.log
├── holdoutCombined.raw
├── holdoutCombinedRaw.log
├── holdoutID.txt
├── merged_allChromosomes.bed
├── merged_allChromosomes.bim
├── merged_allChromosomes.fam
├── merged_allChromosomes.log
├── merged_allChromosomes.snplist
├── combinedID.txt
├── testCombined.raw
├── testCombinedRaw.log
├── testID.txt
├── trainingCombined.raw
├── trainingCombinedRaw.log
├── trainingID.txt
├── pheno.config
├── epiFiles
│   ├── preSummaryFiles
│   │   ├── trainingEpi.epi.cc
│   │   ├── trainingEpi.epi.cc.1
│   │   ├── trainingEpi.epi.cc.10
│   │   ├── trainingEpi.epi.cc.11
│   │   ├── trainingEpi.epi.cc.12
│   │   ├── trainingEpi.epi.cc.13
│   │   ├── trainingEpi.epi.cc.14
│   │   ├── trainingEpi.epi.cc.15
│   │   ├── trainingEpi.epi.cc.16
│   │   ├── trainingEpi.epi.cc.17
│   │   ├── trainingEpi.epi.cc.18
│   │   ├── trainingEpi.epi.cc.2
│   │   ├── trainingEpi.epi.cc.3
│   │   ├── trainingEpi.epi.cc.4
│   │   ├── trainingEpi.epi.cc.5
│   │   ├── trainingEpi.epi.cc.6
│   │   ├── trainingEpi.epi.cc.7
│   │   ├── trainingEpi.epi.cc.8
│   │   ├── trainingEpi.epi.cc.9
│   │   ├── trainingEpi.epi.cc.summary.1
│   │   ├── trainingEpi.epi.cc.summary.10
│   │   ├── trainingEpi.epi.cc.summary.11
│   │   ├── trainingEpi.epi.cc.summary.12
│   │   ├── trainingEpi.epi.cc.summary.13
│   │   ├── trainingEpi.epi.cc.summary.14
│   │   ├── trainingEpi.epi.cc.summary.15
│   │   ├── trainingEpi.epi.cc.summary.16
│   │   ├── trainingEpi.epi.cc.summary.17
│   │   ├── trainingEpi.epi.cc.summary.18
│   │   ├── trainingEpi.epi.cc.summary.2
│   │   ├── trainingEpi.epi.cc.summary.3
│   │   ├── trainingEpi.epi.cc.summary.4
│   │   ├── trainingEpi.epi.cc.summary.5
│   │   ├── trainingEpi.epi.cc.summary.6
│   │   ├── trainingEpi.epi.cc.summary.7
│   │   ├── trainingEpi.epi.cc.summary.8
│   │   ├── trainingEpi.epi.cc.summary.9
│   │   └── trainingEpi.log
│   ├── trainingCombinedEpi.epi.cc.summary
│   ├── trainingCombinedEpi.epi.cc.summary.filtered
│   └── trainingCombinedEpi.log
├── figures
│   ├── AUC_metrics_table_{prs1/prs2..}.{nfeatures}.holdout.{mixed/protect/risk}.csv
│   ├── AUC_metrics_table_{prs1/prs2..}.{nfeatures}.validation.{mixed/protect/risk}.csv
│   ├── {prs1/prs2..}.{nfeatures}.holdout.{mixed/protect/risk}.AUC.png
│   ├── {prs1/prs2..}.{nfeatures}.validation.{mixed/protect/risk}.AUC.png
│   ├── {prs1/prs2..}.{nfeatures}.holdout.{mixed/protect/risk}.boxplot.png
│   ├── {prs1/prs2..}.{nfeatures}.validation.{mixed/protect/risk}.boxplot.png
│   ├── {prs1/prs2..}.{nfeatures}.holdout.{mixed/protect/risk}.densityPlot.png
│   ├── {prs1/prs2..}.{nfeatures}.validation.{mixed/protect/risk}.densityPlot.png
│   ├── {prs1/prs2..}.{nfeatures}.holdout.{mixed/protect/risk}.prevalencePlot.png
│   ├── {prs1/prs2..}.{nfeatures}.validation.{mixed/protect/risk}.prevalencePlot.png
│   ├── importantFeatureZscores.{batchIteration}.cardio.png
│   ├── importantFeatureZscores.{batchIteration}.epi.png
│   ├── importantFeatureZscores.{batchIteration}.main.png
│   ├── shap_summary_plot.{batchIteration}.cardio.png
│   ├── shap_summary_plot.{batchIteration}.epi.png
│   └── shap_summary_plot.{batchIteration}.main.png
├── models
│   ├── imp_mean_{epi/main}_{batchIteration}.pkl
│   ├── sklearnGradBoostHistClassifier_cardioMetabolic_{env1/env2..}.pkl
│   ├── sklearnGradBoostHistClassifier_{epi/main}_{batchIteration}.pkll
│   └── sklearnNaiveBayes_{epi/main}_{batchIteration}.pkl
└─── scores
    ├── cardio.{nFeaturesInPRS}.holdout.{mixed/protect/risk}.prs.csv
    ├── cardio.{nFeaturesInPRS}.validation.{mixed/protect/risk}.prs.csv
    ├── epi.{nFeaturesInPRS}.holdout.{mixed/protect/risk}.prs.csv
    ├── epi.{nFeaturesInPRS}.validation.{mixed/protect/risk}.prs.csv
    ├── epi+main.{nFeaturesInPRS}.holdout.{mixed/protect/risk}.prs.csv
    ├── epi+main.{nFeaturesInPRS}.validation.{mixed/protect/risk}.prs.csv
    ├── main.{nFeaturesInPRS}.holdout.{mixed/protect/risk}.prs.csv
    ├── main.{nFeaturesInPRS}.validation.{mixed/protect/risk}.prs.csv
    ├── all.{nFeaturesInPRS}.holdout.{mixed/protect/risk}.prs.csv
    ├── all.{nFeaturesInPRS}.validation.{mixed/protect/risk}.prs.csv
    ├── cardioMetabolicimportantFeaturesPostShap.csv
    ├── cardioMetabolicModelScores.csv
    ├── featureScores.csv
    ├── importantFeaturesForAssociationAnalysis.csv
    ├── importantFeaturesPostShap.csv
    └── sklearnModelScoresSections.csv


6 directories


```


# WORKFLOW OVERVIEW: 

## G, GxG, and GxGxE Analysis Overview using T2D data as an example

![PRS Pipeline Workflow](READMEfigures/simplifiedWorkflow.png)


# RUN TEST WORKFLOW

## clone github repository to local machine

## INSTALL [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install, "Install instructions")

## Test Workflow overview

## creates and activates a conda environment ukb_env

## creates folders and results for type2Diabetes_test using steps:

- create simulation data in testData directory:
  + genotyped data
  + environmental data
  + hla data
  + covariate data
  
- use simulated data for all steps to calculating PRS for G, GxG, G+(GxG) and GxGxE datasets

- results will be found in the prsInteractive/results/type2Diabetes_test



``` 
cd prsInteractive/workflows/

bash run_workflow_test.sh

```



# RUN HPC WORKLOW 


## setup environment variables and create and activate conda environment


- create hla, environmental, and covariate data in prsInteractive/results/:

+ environmental data
+ hla data
+ covariate data
+ pheno.config file

- create pheno specific directories in results/pheno/

+ scores/
+ figures/
+ models/
+ epiFiles/

* pheno = phenotype spelled in camel font and no spaces (i.e. type2Diabetes)
* icd10 code = substring present in the UK Biobank data (i.e. E11)
* pheno substring = will be exact spelling found in UKB data to check for if icd10 not present. this will have spaces so will need to wrap in " " (i.e. "type 2 diabetes")
* n = number of cores to pass to the epistatic analysis, with 40 cores being the norm and will take approximately 48 hours (i.e. 20)


```bash

$ cd /path/to/directory/prsInteractive


$ bash ../envSetUp.sh $pheno $icd10 "${phenoStr}" $n (# of cores on local machine)


```

## run the workflow in order

### 1) cleans variants, starts fast-epistasis analysis using plink and begins batch feature reduction step for main (single SNPs)
  
  ```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_data_cleaning_workflow_submit.sh {pheno} {icd10 code}  {"sub string"} {n}
  
  ```
  

### 2) when fast-epistasis is complete, run the batch feature reduction for epi features

  
  ```bash 
  
$ cd path/to/prsInteractive/hpc
  
$ sbatch run_model_batches_submit.sh pheno "epi"
  
  ````

### 3) when feature reduction is complete for main and epi, run the GxGxE interaction analysis

#### envStr = user decision which will be saved as ENV_TYPE and used in file names (i.e. cardioMetabolic was used in thesis)

```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_gene_environment_feature_discovery_submit.sh pheno "envStr"

````

### 4) create combined EnvGeno matrix to be used downstream 

```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_create_gene_env_data_submit.sh pheno 

````

### 4) create combined EnvGeno matrix to be used downstream 

```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_create_gene_env_data_submit.sh pheno 

````

### 5) run final association analysis with reduced features 

```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_glmNetFinalModel.sh pheno 

````

### 6) calculate PRS for all models trained on individual and combined features sets 

```bash 

$ cd path/to/prsInteractive/hpc

$ sbatch run_create_gene_env_data_submit.sh pheno 

````


  
# Running analysis with WDL workflow
  
########## WORKING IN PROGRESS ########

### These instructions use the cromwell backend to compile .wdl workflow
###
### cromwell specific features include a cromwell.config file and command line instruction to start run
###

### Input

All inputs are required for from the root directory and annotated in the pipelineInputs.json file
Phenotype specific inputs are updated in first step
###

### Download cromwell 

Instructions for download can be found here: [cromwell download] (https://cromwell.readthedocs.io/en/latest/tutorials/FiveMinuteIntro/)



#### Step 1: Setup Environment needed to run workflow

#### inputs:
  * config/default.config
  * data/
    - covar.txt
    - hla_participant.csv
    - participant.csv
    - participant_environment.csv
    - ukb_hla_V2.txt
    - withdrawals.csv

#### this step creates:
* .env 
* pipelineInputs.json
* directories for the phenotype
  - results/<pheno>
* results/<pheno>/summary.txt #file of inputs used
* results/<pheno>/pheno_config.sh file

```
cd /prsInteractive

tar -czf scripts.tar.gz scripts

./envSetUp.sh <phenotype> <icd10> <phenotype string used in search> <n cores for epistatic analysis> <platform: (hpc,local,dnanexus) to run analysis>

```
  
  
  



  

    