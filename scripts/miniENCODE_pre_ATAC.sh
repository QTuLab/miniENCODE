#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes ATAC-seq data by aligning with Bowtie2 and calling peaks using MACS2.

# Define common variables
dir=$PWD
dataDir=$dir/data
scriptDir=$dir/scripts
infoDir=$dir/data/info
thread=4

# Step 1: ATAC-seq Alignment with Bowtie2
echo "Mapping by Bowtie2"
python $scriptDir/miniENCODE_pre_ATAC_alignment.py -i $infoDir/id_atac.txt -o $dataDir -t $thread 

# Step 2: Calling Peaks
echo "Calling Peaks"
python $scriptDir/miniENCODE_pre_ATAC_callpeak.py -s $infoDir/namelist_atac.csv -o $dataDir -t $thread 
