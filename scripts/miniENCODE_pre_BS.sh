#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes BS-seq data by mapping with Bismark and merging the resulting data.

# Define common variables
dir=$PWD
dataDir=$dir/data
scriptDir=$dir/scripts
infoDir=$dir/data/info
thread=4

# Step 1: Mapping by Bismark
echo "Mapping by Bismark"
time nice python $scriptDir/miniENCODE_pre_BS_alignment.py -i $infoDir/id_bs.txt -o $dataDir -t $thread 

# Step 2: Merge Data
echo "Merging Data"
time nice python $scriptDir/miniENCODE_pre_BS_rename_bw.py -s $infoDir/namelist_bs.csv -o $dataDir -t $thread 

rm -rf ${scriptDir}/__pycache__
