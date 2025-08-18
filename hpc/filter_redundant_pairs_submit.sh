#!/bin/bash

#SBATCH --job-name=model_batches
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=keri@multerer.com
#

module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env


EPI_FILE=$1
SCRIPTS_DIR=$2
export EPI_FILE=$EPI_FILE

python "$SCRIPTS_DIR/filter_redundant_epi_pairs.py"

echo "Python script filtering redundant epi features finished with exit code: $?"