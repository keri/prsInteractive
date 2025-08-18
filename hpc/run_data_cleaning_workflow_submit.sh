#!/bin/bash -x

#
#SBATCH --job-name=clean_create_variant_data
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_clean_create_variant_data.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_clean_create_variant_data.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=70G
#SBATCH --time=7:00:00
#

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   raw variant calls with bed format : ukb{project_number}_c1_b0_v2
#   raw hla data in csv format  : hla_participant.csv',index_col='Participant ID') and file with headers : ukb_hla_v2.txt
#   participant data in csv format with 
platform=${5:-"local"} 
pheno=$1
icd10=$2
phenoStr=$3
N_CORES=$4

#pheno="myocardialInfarction"



module load plink/1.90



##############  SET UP ENV VARIABLES FOR JOB #################

# Check if parent script exists
if [[ ! -f "../envSetUp.sh" ]]; then
    echo "Error: envSetUp.sh script not found at ../envSetUp.sh"
    exit 1
fi

# Check if envSetUp has already been run
if [[ ! -d "../results/$pheno" ]]; then
    echo "Error: envSetUp.sh script not run"
    # Run parent script with all arguments
    echo "Executing: ../envSetUp.sh"
    cd ../
    bash envSetUp.sh $pheno $icd10 "${phenoStr}" $N_CORES
    cd hpc/
fi

# Capture exit status
exit_status=$?
if [[ $exit_status -ne 0 ]]; then
    echo "envSetUp.sh failed with exit status: $exit_status"
    exit $exit_status
fi

echo "envSetUp.sh completed successfully"

module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env


# Source config
source ../env.config # because you're in prsInteractive/hpc


#check that a results folder for phenotype exists
if [ ! -d "${RESULTS_PATH}/$pheno" ]; then
    echo "Folder '${RESULTS_PATH}/$pheno' does not exist..."
    echo "run envSetUp.sh <pheno> <icd10> <phenoStr> <n cores to use in epistatic interaction analysis>"
    exit 1
    
else
    echo "sourcing $pheno env variables."
    #source pheno specific environment variables
    source "${RESULTS_PATH}/$pheno/pheno.config"
fi

export PHENO_PATH=$PHENO_PATH
export PHENO=$PHENO
export PHENO_STR=$PHENO_STR
export ICD10=$ICD10
export EPI_PATH=$EPI_PATH
export DATA_PATH=$DATA_PATH
export RESULTS_PATH=$RESULTS_PATH
export N_CORES=$N_CORES
export PRS_INTERACTIVE_HOME=$PRS_INTERACTIVE_HOME
export WITHDRAWAL_PATH=$WITHDRAWAL_PATH


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: $PHENO_STR"

#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_PATH}" ]; then
    echo "Folder '${PHENO_PATH}' does not exist. Creating it..."
    mkdir "${PHENO_PATH}" 
    
else
    echo "Folder '${PHENO_PATH}' already exists."	
fi




###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############

# create phenotype data and train test split IDs
python "${SCRIPTS_DIR}/create_pheno_train_test_split.py"

# create hla (and environmental data files?)
python "${SCRIPTS_DIR}/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
bash "${SCRIPTS_DIR}/plink_clean_variant_calls.sh"

#merge separate chromosome files into one
bash "${SCRIPTS_DIR}/merge_chromosomes.sh"

sbatch multiprocessing_fast_epistasis_submit.sh

sbatch run_model_batches_submit.sh $PHENO "main"










