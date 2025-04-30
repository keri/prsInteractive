#python script keri Multerer October 2022
#script to clean .raw column headings from rsID_A/C/T/G to rsID

import pandas as pd
import os
import argparse

def get_columns(inputfilepath):
    """input : chromosome#,mainfilepath
    output: text file with space separated column headings cleaned"""

    df = pd.read_csv(inputfilepath, sep=" ",nrows=1)

    return (df.columns.to_list())


def create_column_file(outputfilepath):
    """input : filepath = absolute filepath to .map
               raw = [FID,IID,PAT,MAT,SEX,PHENOTYPE,snpId1_(effect allele),snpId2__(effect allele).....snpIdk_(effect allele)]
               cells = [0,1,2,NaN]
               .map = columns=[0,1,2,3] with 1 being the rsIDs
        ouput : columns = same as input except rsID_MAF to rsID"""

    for file in ["merged_allChromosomes.raw","holdoutCombined.raw"]:
        inputFile = f"{outputfilepath}/{file}"
        snp_columns = get_columns(inputFile)

        cleaned_columns = [col.split('_')[:-1] for col in snp_columns]

        columns_file = pd.DataFrame(columns=cleaned_columns)
        columns_file.to_csv(outputfilepath+f"/{file}_columns.txt",sep=" ",index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="cleaning columns ...")
    parser.add_argument("--data_folder", help="Path to the input data folder")
    parser.add_argument("--pheno_folder", help="Path to the input pheno data folder")
    parser.add_argument("--pheno", help="Phenotype to analyze")
    parser.add_argument("--pheno_str", help="Phenotype string used in Non-cancer illness code, self-reported field")
    parser.add_argument("--icd_code", help="ICD 10 code of phenotype")
    
    
    args = parser.parse_args()
    
    # Prefer command-line input if provided; fallback to env var
    data_path = args.data_folder or os.environ.get("DATA_PATH")
    print(f"[PYTHON] Reading from: {data_path}")
    
    pheno_path = args.pheno_folder or os.environ.get("PHENO_PATH")
    print(f"[PYTHON] Reading from: {pheno_path}")

    mainfilepath = f'/nfs/scratch/multerke/ukbiobank/{pheno}'
    #mainfilepath = '/nfs/scratch/multerke/ukbiobank/bulkFiles/cleaned'

    create_column_file(mainfilepath,pheno)
    