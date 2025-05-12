#Steps to run pipeline on the hpc


##Run the workflow to clean and get the fast-epistasis analysis going

###Input files must be present:
   in prsInteractive/data/
    1. covar.txt 
    2. hla_participant.csv
    3. participant.csv
    4. participant_environment.csv
    5. ukb_hla_v2.txt
    6. variant_calls/
        genotyped data in .bed/.bim/.fam
    7. withdrawals.csv
    8. withdrawalsID.txt
    
###Output files in pheno folder:
    ####/nfs/scratch/projects/ukbiobank/prsInteractive/results/{pheno}/
    1. new phenotype folder created for user input pheno argument (for run_data_cleaning_workflow_submit.sh)
    2. combinedID.txt
    3. holdoutID.txt
    4. testID.txt
    5. trainingID.txt
    6. merged_allChromosomes.bed/.bim/.fam
    7. merged_allChromosomes.raw
    8. holdoutCombined.bed/.bim/.fam
    9. holdoutCombined.raw
    10. merged_allChromosomes.snplist
    11. pheno.txt
    12. epiFiles/
        preSummaryFiles/
    ####/nfs/scratch/projects/ukbiobank/prsInteractive/results/
    1. participant_environment.csv
    2. participant_hla.csv
    3. covar.txt

#########################  STEPS IN ANALYSIS  #############################

###1) check to see if a conda environment has been created in /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
 - if ukb_env is not present, create it with this command:
    # Create environment from file
    # Option 1: Create ukb_env environment set in the environment.yml file
    # this will create a conda "ukb_env" folder with all of the dependencies in the nfs/scratch/projects/ukbiobank/prsInteractive/ directory
    # run this command
    
    username@raapoi-login:/nfs/scratch/projects/ukbiobank/prsInteractive$ conda env create -f environment.yml
    
###2) with conda env "ukb_env" present run the workflow:
    #run command lines:
    
    $ cd nfs/scratch/projects/ukbiobank/prsInteractive/hpc
    
    $ sbatch run_data_cleaning_workflow_submit.sh {pheno name you create} {icd10 code} substring present in the UK Biobank data {"sub string"} {n cores for epistatic interaction run, normally used 40} 
    
        i.e. sbatch run_data_cleaning_workflow_submit.sh myocardialInfarction I21 "myocardial infarction" 40
    
    #has the following scripts
    # create phenotype data and train test split IDs
    python "${SCRIPTS_DIR}/create_pheno_train_test_split.py"
    
    # create hla (and environmental data files?)
    python "${SCRIPTS_DIR}/clean_environment_hla_covar_data.py"
    
    # Run the variant call cleaning
    bash "${SCRIPTS_DIR}/plink_clean_variant_calls.sh"
    
    #merge separate chromosome files into one
    bash "${SCRIPTS_DIR}/merge_chromosomes.sh"
    
    sbatch multiprocessing_fast_epistasis_submit.sh
    

    
    
