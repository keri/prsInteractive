#!/bin/bash

#
#SBATCH --job-name=clean_create_variant_data
#SBATCH -o  /nfs/scratch/projects/ukbiobank/ashley/err_out/%A_clean_create_variant_data.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/ashley/err_out/%A_clean_create_variant_data.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=70G
#SBATCH --time=2:00:00
#

################## USE DATA FROM UK BIOBANK ############
# FILES THAT MUST BE PRESENT:
#   raw variant calls with bed format : ukb{project_number}_c1_b0_v2
#   raw hla data in csv format  : hla_participant.csv',index_col='Participant ID') and file with headers : ukb_hla_v2.txt
#   participant data in csv format with 

pheno=$1
icd10=$2
phenoStr=$3
n=$4

#pheno="myocardialInfarction"
#icd10="I21"
#phenoStr="myocardial infarction"

module load Miniconda3/4.9.2
eval "$(conda shell.bash hook)"
conda activate /nfs/home/multerke/.conda/envs/py39

module load plink/1.90

# Set base directories
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$WORKFLOW_DIR")"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
PHENO_DIR="$RESULTS_DIR/$pheno"
HPC_DIR="$PROJECT_ROOT/hpc"


# Set environment variable
export DATA_PATH="$DATA_DIR"
export RESULTS_PATH="$RESULTS_DIR"
export PHENO_PATH="$PHENO_DIR"
export PHENO="$pheno"
export PHENO_STR="$phenoStr"
export ICD="$icd10"
export N=$n


echo "[WORKFLOW] DATA_PATH is set to: $DATA_PATH"
echo "[WORKFLOW] RESULTS_PATH is set to: $RESULTS_PATH"
echo "[WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "[WORKFLOW] Scripts directory: $SCRIPTS_DIR"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "ICD 10 BEING ANALYZED ...: $ICD"
echo "PHENOTYPE STRING TO FILTER FOR IF ICD CODE NOT PRESENT ...: $PHENO_STR"




#make a folder inside root data folder for each phenotype

if [ ! -d "${PHENO_DIR}" ]; then
    echo "Folder '${PHENO_DIR}' does not exist. Creating it..."
    mkdir "${PHENO_DIR}" 

    
else
    echo "Folder '${PHENO_DIR}' already exists."	
fi


# create phenotype data and train test split IDs
python "$SCRIPTS_DIR/create_pheno_train_test_split.py"

# create hla (and environmental data files?)
python "$SCRIPTS_DIR/clean_environment_hla_covar_data.py"

# Run the variant call cleaning
bash "$SCRIPTS_DIR/plink_clean_variant_calls.sh"

#merge separate chromosome files into one
bash "$SCRIPTS_DIR/merge_chromosomes.sh"

sbatch "$HPC_DIR/multiprocessing_fast_epistasis_submit.sh"










