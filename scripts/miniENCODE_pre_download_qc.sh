#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Download and perform quality control on sequencing data.

# Define directory paths
dir=$PWD
dataDir="${dir}/data/"
infoDir="${dir}/data/info/"
scriptDir="${dir}/scripts/"

# Set the number of threads
thread=4

# Message indicating the start of the download and quality control process
echo 'Download and Quality Control'

# Execute the Python script for downloading and quality control
python "${scriptDir}/miniENCODE_pre_download_qc.py" -i "${infoDir}/id_download.txt" -o "${dataDir}" -t "${thread}"


