#!/bin/bash


#
#SBATCH --job-name=linkage_disequilibrium
#SBATCH -o  /nfs/scratch/projects/ukbiobank/err_out/%A_LD.out
#SBATCH -e /nfs/scratch/projects/ukbiobank/err_out/%A_LD.err
#SBATCH --partition=quicktest
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=00:45:00
#


pheno=$1

#pheno="myocardialInfarction"




module load Miniconda3/4.9.2
source $(conda info --base)/etc/profile.d/conda.sh 
conda activate /nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env
#export PATH="/nfs/scratch/projects/ukbiobank/prsInteractive/ukb_env/bin:$PATH"

module load plink/1.90


echo "PHENO is set to : $PHENO"



# Source config with error handling
if [ ! -f "../env.config" ]; then
	echo "ERROR: ../env.config not found!"
	echo "Current directory: $(pwd)"
	echo "Looking for: $(realpath ../env.config 2>/dev/null || echo '../env.config')"
	exit 1
	
else
	source ../env.config
fi

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

echo "[DIR] scripts directory : $SCRIPTS_DIR"

export PHENO_PATH=$PHENO_PATH

python "$SCRIPTS_DIR/helper/create_LD_snpList.py"

wait 

plink --bfile "$PHENO_PATH/merged_allChromosomes" --extract "$PHENO_PATH/finalModelLDSnps.txt" --indep-pairwise 100kb 1 .6 --r2 --show-tags all --out "$PHENO_PATH/finalModel"

#overwites the reducedModelFeatureScores.csv with pruned features in LD and reversed
python "$SCRIPTS_DIR/filter_features_LD.py"