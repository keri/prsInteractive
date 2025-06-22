#Steps to run pipeline on the hpc


##Run the workflow to clean and get the fast-epistasis analysis going

###Input files must be present:
   in prsInteractive/data/
    1. covar.csv
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
    3. covar.csv

#########################  STEPS IN ANALYSIS  #############################

###1) check to see if a conda environment has been created in /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
 - if ukb_env is not present, create it with this command:
    # Create environment from file
    # Option 1: Create ukb_env environment set in the environment.yml file
    # this will create a conda "ukb_env" folder with all of the dependencies in the nfs/scratch/projects/ukbiobank/prsInteractive/ directory
    # run these commands preceded with "$"
    
    $ cd /nfs/scratch/projects/ukbiobank/prsInteractive
    $ module load Miniconda3/23.9.0-0
    
    Your command line prompt should look like this:
    username@raapoi-login:/nfs/scratch/projects/ukbiobank/prsInteractive$ 
    
    $ conda env create --prefix ./ukb_env -f environment.yml
    
    
###2) with conda env "ukb_env" present run the workflow:

    pheno = phenotype spelled in camel font and no spaces (i.e. type2Diabetes, myocardialInfarction)
    icd10 code = substring present in the UK Biobank data
    pheno substring = will be exact spelling found in UKB data to check for if icd10 not present. this will have spaces so will need to wrap in " "
    n = number of cores to pass to the epistatic analysis, with 40 cores being the norm and will take approximately 48 hours
    
    #run command lines:
    
    $ cd nfs/scratch/projects/ukbiobank/prsInteractive/hpc
    
    $ sbatch run_data_cleaning_workflow_submit.sh {pheno} {icd10 code}  {"sub string"} {n} 
    
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
    
####3) After epistatic analysis is complete, run the batch models with 
    from the hpc/directory:
    pheno = name of folder created and entered in run_data_cleaning_workflow_submit.sh
    data_type = main (if running single SNPs) or epi (if running with epi-pairs created from epistatic analysis
    
    $ sbatch run_model_model_batches_submit.sh {pheno} {data_type}
    

    
    
    