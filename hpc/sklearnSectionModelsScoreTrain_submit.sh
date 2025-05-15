#!/bin/bash

#!/bin/bash

#
#SBATCH --job-name=train_score_section_models
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=50
#SBATCH --mem=700G
#SBATCH --time=05:00:00
#

module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"

# Source config
#source ../config.sh  # because you're in prsInteractive/hpc


echo "[HPC WORKFLOW] PHENO_PATH is set to: $PHENO_PATH"
echo "starting iteration ...: $START"
echo "ending iteration ...: $END"
echo "PHENOTYPE BEING ANALYZED ...: $PHENO"
echo "[HPC WORKFLOW] SCRIPTS_PATH is set: $SCRIPTS_DIR"
echo "[HPC WORKFLOW] TRAINING_PATH is set: $TRAINING_PATH"
echo "[HPC WORKFLOW] TEST_PATH is set: $TEST_PATH"


###########  CREATE A PHENOTYPE FOLDER TO COLLECT RESULTS IF NOT PRESENT ############


python "${SCRIPTS_DIR}/sklearnSectionModelsScoreTrain.py"











