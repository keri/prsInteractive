#!/bin/bash

#!/bin/bash

#
#SBATCH --job-name=train_score_main__models
#SBATCH -o  /nfs/scratch/multerke/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/multerke/ukbiobank/err_out/%A.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=50
#SBATCH --mem=700G
#SBATCH --time=4:00:00
#

module load Miniconda3/23.9.0-0
eval "$(conda shell.bash hook)"
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"

# Source config
#source ../config.sh  # because you're in prsInteractive/hpc


echo "[HPC WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "starting iteration ...: $START"
echo "ending iteration ...: $STOP"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_PATH"
echo "[HPC WORKFLOW] TRAINING_PATH is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST_PATH is set: $TEST_PATH"


###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


python "${SCRIPTS_PATH}/sklearnSectionModelsScoreTrain.py"











