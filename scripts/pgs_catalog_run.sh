#!/bin/bash

phenoPath=$1


conda activate pgs_tools

#type2Diabetes prs
#nextflow run pgscatalog/pgsc_calc -profile docker --input "${phenoPath}/pgs_catalog_samplesheet.csv" --pgs_id PGS000020,PGS001357 --target_build GRCh37 --output "${phenoPath}/pgps_catalog_scores/" --min_overlap 0

#celiacDisease prs
nextflow run pgscatalog/pgsc_calc -profile docker --input "${phenoPath}/pgs_catalog_samplesheet.csv" --pgs_id  --target_build GRCh37 --output "${phenoPath}/pgps_catalog_scores/" --min_overlap 0