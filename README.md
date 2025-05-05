# OVERVIEW

The prs Insteractive pipeline was developed using genotyped, imputed HLA, and environmental data for UK Biobank participants with a European background. The pipeline is comprised of workflows for the inclusion of inclusion of gene (G), gene-gene (GxGxE), and gene-environment (GxE) interaction weights into polygenic risk (PRS) calculations for complex traits. The pipeline was developed and validated for type 2 diabetes (T2D) and celiac disease (CD) but can be applied to any trait with a ICD10 code and/or substring filter present in the Non-cancer illness code [instances 0-2] fields. Specific G, GxG, and GxGxE cohorts are further analysed to identify underlying features important to each cohorts to identify different molecular pathways driving risk within the cohorts.


#INPUT FILES NEEDED 

## European ancestry filters used in development dataset
### Fields:
 Used in genetic principal components = Yes
 Outliers for heterozygosity or missing rate = NULL
 Sex chromosome aneuploidy = NULL
 Genetic kinship or other participants = NOT ANY OF Ten or more third-degree relatives identified
 Genetic ethnic groups = IS NOT NULL
 Genetic principal components | Array 1 = IS BETWEEN -20 - 40
 Genetic principal components | Array 2 = IS BETWEEN -25 - 10

## participant.csv 
File downloaded from cohort created using cohort browser on DNA nexus platform.

### minimum columns needed : 
 ["Participant ID",
  "Diagnoses - main ICD10",
  "Non-cancer illness code, self-reported | Instance 0",
    "Non-cancer illness code, self-reported | Instance 1",
    "Non-cancer illness code, self-reported | Instance 2",
    "Non-cancer illness code, self-reported | Instance 3",
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
   "Age at recruitment"]
    

## raw variant calls in bed format

## hla_participant.csv

### column headings = ukb_hla_v2.txt 
### imputed data = hla_participant.csv
File downloaded from cohort created using cohort browser on DNA nexus platform with imputation values range from 0 - 2

## participant_environment.csv
 File downloaded from cohort created using cohort browser on DNA nexus platform with clinical markers for blood chemistry, blood counts, and cardiometabolic features listed in Supplemental Table S13 of thesis
 
 Workflow can use as input any clinical marker which includes blood counts, blood chemistry, and lifestyle data available at initial screening and follow up visits for all participants.
    

## withdrawals.csv
 A list of eid's provided by UK Biobank of people who have opted out of research. File consists of one column with no heading or index.


# G, GxG, and GxGxE Analysis Overview using T2D data as an example

[insert workflow image : combinedGWASWorkflow.png]

# Important Underlying Feature Analysis using T2D data as an example

[insert important feature workflow image : importantFeatureWorkflowSHAP.png]

# File Structure needed for running analysis



