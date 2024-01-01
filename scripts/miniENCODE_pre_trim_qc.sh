#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Perform trimming and quality control on various sequencing data types.

# Define common variables
dir=$PWD
dataDir="${dir}/data/"
infoDir="${dir}/data/info/"
scriptDir="${dir}/scripts/"
thread=4

echo "Trimming and Quality Control"

# Loop through different ID files for RNA-seq, ATAC-seq, ChIP-seq, and BS-seq
for idtxt in id_rna.txt id_atac.txt id_chip.txt id_bs.txt
do
    echo "Processing ${idtxt}"
    python "${scriptDir}/miniENCODE_pre_trim_qc.py" -i "${infoDir}/${idtxt}" -o "${dataDir}" -t "${thread}"
done
