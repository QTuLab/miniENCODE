#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Process and finalize the enhancer-like signature (ELS) data.

# Define the data directory
dir=$PWD
dataDir="${dir}/data/ELS/"

# Create a temporary directory for intermediate files
mkdir -p "${dataDir}/tmp"
cd "${dataDir}/tmp" || exit
cp "${dataDir}/atac_peaks.bed" .

# Step 1: Identify ELS based on peak cutoffs
echo "Step 1 ..." >> els.log
minP=4942

cp atac_peaks.bed 1
for j in `seq 2 1 $minP`
do
    cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
    awk -v a="$j" '{if($9>=a) {print $0}}' 1 > 2
    bedtools merge -d 1 -c 5 -o min -i 2 | awk '{if ($3-$2 >= 50) print $0}' > "els_cutoff.$j.bed"
    mv -f 2 1
    num=$(wc -l "els_cutoff.$j.bed" | awk '{print $1}')
    echo -e "\t$cutoff $num" >> els.log
done

# Step 2: Refine ELS by merging overlapping regions
echo "Step 2 ..." >> els.log
awk '{if ($3-$2 <= 350) print $0}' els_cutoff.2.bed > peaks
for j in `seq 3 1 $minP`
do
    cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
    bedtools intersect -v -a "els_cutoff.$j.bed" -b peaks > tmp
    awk '{if ($3-$2 <= 500) print $0}' tmp >> peaks
    num=$(wc -l tmp | awk '{print $1}')
    echo -e "\t$cutoff $num" >> els.log
done

mv -f peaks at_cels.bed
cp at_cels.bed "${dataDir}/rDAR_process.bed"

# Move to the data directory
cd "${dataDir}" || exit

# Finalize the ELS data
awk 'BEGIN {rank=0; before=0; running=1}{if ($2 != before) rank = running; print $1 "\t" $2 "\t" $3 "\tpeak_" rank "\t" $4; before=$2; running += 1}' "${dataDir}/rDAR_process.bed" > "${dataDir}/rDAR_final.bed"

# Clean up the temporary directory
rm -rf "${dataDir}/tmp"
