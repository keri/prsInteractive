#!/bin/bash -x

#
#SBATCH --job-name=train_score_section_models
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=60
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#

#module load Miniconda3/4.9.2
#source $(conda info --base)/etc/profile.d/conda.sh 
#conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"

# Source config
#source ../env.config  # because you're in prsInteractive/hpc


echo "[HPC WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "starting iteration ...: $START"
echo "ending iteration ...: $END"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_DIR"
echo "[HPC WORKFLOW] TRAINING_PATH is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST_PATH is set: $TEST_PATH"
echo "[EPI_FILE EPI_PATH is set: $EPI_FILE"
echo "DATA_TYPE is set: $DATA_TYPE"


export START=$START
export END=$END
export DATA_TYPE=$DATA_TYPE
export PHENO=$PHENO 
export PHENO_PATH=$PHENO_PATH
export TRAINING_PATH=$TRAINING_PATH
export TEST_PATH=$TEST_PATH
export EPI_FILE=$EPI_FILE



python "${SCRIPTS_DIR}/sklearnSectionModelsScoreTrain.py"

#conda deactivate









