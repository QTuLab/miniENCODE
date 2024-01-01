#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Process RNA-seq data including mapping with HISAT2, merging, and quantification using StringTie.

# Define common variables
dir=$PWD
dataDir="${dir}/data"
infoDir="${dir}/data/info"
scriptDir="${dir}/scripts"
thread=4

# Step 1: RNA Alignment with HISAT2
echo "Mapping by HISAT2"
python "${scriptDir}/miniENCODE_pre_RNA_alignment.py" -i "${infoDir}/id_rna.txt" -o "${dataDir}" -t "${thread}"

# Step 2: Merge Data
echo "Merging Data"
python "${scriptDir}/miniENCODE_pre_RNA_merge.py" -s "${infoDir}/namelist_rna.csv" -o "${dataDir}" -t "${thread}"

# Step 3: Quantify by StringTie
echo "Quantify by StringTie"
python "${scriptDir}/miniENCODE_pre_RNA_quantify.py" -s "${infoDir}/namelist_rna.csv" -o "${dataDir}" -t "${thread}"
