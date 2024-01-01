#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes ChIP-seq data by aligning with Bowtie2 and calling peaks.

# Define common variables
dir=$PWD
dataDir=$dir/data
scriptDir=$dir/scripts
infoDir=$dir/data/info
thread=4

# Step 1: Mapping by Bowtie2
echo "Mapping by Bowtie2"
python $scriptDir/miniENCODE_pre_ChIP_alignment.py -i $infoDir/id_chip.txt -o $dataDir -t $thread 

# Step 2: Calling Peaks
echo "Calling Peaks"
python $scriptDir/miniENCODE_pre_ChIP_callpeak.py -s $infoDir/namelist_chip.csv \
                  -c $infoDir/namelist_chip_callpeak.csv -o $dataDir -t $thread 
