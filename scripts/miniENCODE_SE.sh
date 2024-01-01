#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Identify Super Enhancers using ROSE.

# Set the path to ROSE
dir=$PWD
PATHTO="/usr/local/share/ROSE-master"
PYTHONPATH="$PATHTO/lib"
export PYTHONPATH
export PATH="$PATH:$PATHTO/bin"
ROSE="$PATHTO/bin/ROSE_main.py"

# Set the number of threads to use
threads=2

# Define the working directories
inputDir="${dir}/data/alignment/"
outDir="${dir}/data/SE/"

# Change to the output directory
cd "$outDir" || exit
chmod +x "$ROSE"
mkdir -p "${outDir}/tmp"

# Define the paths to BAM files and narrowPeak file for the current sample
sp="Heart"
bam="${inputDir}/ch_H3K27ac_heart_2020nat_chr1.final.bam"
inputbam="${inputDir}/ch_input_heart_2020nat_chr1.final.bam"
narrowPeakFile="${dir}/data/peak/ch_H3K27ac_heart_2020nat_chr1_peaks.narrowPeak"
refseq="${dir}/data/reference/GRCz11_ucsc.refseq"

# Reheader BAM files to replace chromosome names
samtools view -@ "$threads" -H "$bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "$bam" > "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam"
samtools index -@ "$threads" "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam"

samtools view -@ "$threads" -H "$inputbam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
    | samtools reheader - "$inputbam" > "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam"
samtools index -@ "$threads" "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam"

# Modify the narrowPeak file to have "chr" prefixes
sed 's/^/chr/g' "$narrowPeakFile" > "${outDir}/tmp/ch_H3K27ac_${sp}_narrowPeak.bed"

# Run ROSE with the specified parameters
mkdir -p "${outDir}/${sp}"
python "$ROSE" --custom="$refseq" -i "${outDir}/tmp/ch_H3K27ac_${sp}_narrowPeak.bed" \
    -r "${outDir}/tmp/ch_H3K27ac_${sp}_chr.bam" \
    -c "${outDir}/tmp/ch_H3K27ac_${sp}_input_chr.bam" \
    -o "${outDir}/${sp}/" \
    -s 12500 -t 2500 2>"${outDir}/${sp}/${sp}.log"

# Remove temporary directory
rm -rf "${outDir}/tmp"
