#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Perform gene regulatory network (GRN) analysis using ANANSE.

# Define reference directory and input/output directories
dir=$PWD
refDir="${dir}/data/reference/"
inputDir="${dir}/data/GRN/input/"
outDir="${dir}/data/GRN/"
reftfbs="${refDir}/GRCz11_TFDB/danRer11.gimme.vertebrate.v5.0.pfm"

# Set the number of threads
thread=1

for sp in Blood
do
    mkdir -p "${outDir}/${sp}"

    # Define input files for ANANSE
    atacbam="${inputDir}/at_${sp}.bam"
    chipbam="${inputDir}/ch_H3K27ac_${sp}.bam"
    elsfile="${inputDir}/els_${sp}.bed"
    tpmfile="${inputDir}/tpm_${sp}.txt"

    # Run ANANSE binding analysis
    ananse binding -A "$atacbam" \
                   -H "$chipbam" \
                   -r "$elsfile" \
                   -g "${refDir}/GRCz11.fa" \
                   -p "$reftfbs" \
                   -o "${outDir}/${sp}" \
                   -n "$thread"

    # Convert ANANSE binding output to TSV format
    ananse view "${outDir}/${sp}/binding.h5" -o "${outDir}/${sp}/binding.tsv"

    # Run ANANSE network analysis
    ananse network -g "${refDir}/GRCz11.fa" \
                   -a "${refDir}/GRCz11.annotation.gene.ensembl.bed" \
                   -e "${tpmfile}" \
                   -o "${outDir}/${sp}/${sp}.network.txt" \
                   -n "$thread" \
                   "${outDir}/${sp}/binding.h5" > "${outDir}/${sp}/log" 2>&1 &
done
