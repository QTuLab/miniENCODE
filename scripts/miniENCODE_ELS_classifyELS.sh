#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Classify and organize enhancer-like signatures (ELS) in ATAC-seq and ChIP-seq data.

# Define the data directory
dir=$PWD
dataDir="${dir}/data/ELS"
refDir="${dir}/data/reference"
scriptDir="${dir}/scripts"
output="${dataDir}/process"

# Set Python environment
pyenv global system

# Define the TSS and BED files
tss="${refDir}/GRCz11_TSS_Basic.bed"
prox="${refDir}/GRCz11_TSS_Basic4K.bed"
rdhs="${dataDir}/rDAR_final.bed"

# Loop through sample names
for sample in Blastula Blood
do
    echo "Processing sample: $sample"
    H3K27ac="${output}/${sample}/H3K27ac"
    cd "${output}/${sample}/ATAC/" || exit

    echo "Splitting ELS into groups..."
    awk '{if ($2 > 1.64) print $0}' maxZ.txt > list
    awk 'FNR==NR {x[$1];next} ($4 in x)' list "$rdhs" > bed
    bedtools intersect -u -a bed -b "$tss" > tss
    bedtools intersect -v -a bed -b "$tss" > a1
    bedtools intersect -u -a a1 -b "$prox" | sort -k1,1 -k2,2n > prox
    bedtools intersect -v -a bed -b "$prox" > distal

    # Calculate center-to-TSS distances
    bedtools closest -d -a prox -b "$tss" > tmp
    cut -f1-4,6-12 tmp > tmp2
    python "${scriptDir}/miniENCODE_ELS_center_distance.py" tmp2 agnostic > new
    awk '{if ($2 >= -200 && $2 <= 200) print $0}' new > center-distance
    awk '{if ($2 < -2000 || $2 > 2000) print $0}' new > far
    awk 'FNR==NR {x[$1];next} ($4 in x)' center-distance prox >> tss
    awk 'FNR==NR {x[$1];next} ($4 in x)' far prox >> distal
    cat center-distance far > new
    awk 'FNR==NR {x[$1];next} !($4 in x)' new prox > tmp
    mv tmp prox

    # Split ELS into different categories
    awk 'FNR==NR {x[$4];next} ($1 in x)' tss "$H3K27ac/maxZ.txt" | \
        awk '{if ($2 > 1.64) print $0}' > pELS
    awk 'FNR==NR {x[$4];next} ($1 in x)' prox "$H3K27ac/maxZ.txt" | \
        awk '{if ($2 > 1.64) print $0}' >> pELS
    awk 'FNR==NR {x[$4];next} ($1 in x)' distal "$H3K27ac/maxZ.txt" | \
        awk '{if ($2 > 1.64) print $0}' > dELS

    # Create accessioned ELS BED files
    awk 'FNR==NR {x[$1];next} ($4 in x)' pELS "$rdhs" | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "pELS"}' > m.bed
    awk 'FNR==NR {x[$1];next} ($4 in x)' dELS "$rdhs" | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "dELS"}' >> m.bed
    
    cat m.bed | sort | uniq > m2.bed
    sort -k1,1 -k2,2n m2.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "_" $4 }' > l.bed
    awk '{$5="1"}1' l.bed > "l2_${sample}.bed"
    cp "${output}/${sample}/ATAC/l.bed" "${dataDir}/${sample}/ELS_${sample}.bed"
done

# Combine all ELS BED files into one aggregated file
cat "${dataDir}"/*/ELS_*.bed | sort -k1,1n -k2,2n -k3,3n | uniq > "${dataDir}/ELS_all_agg_demo.bed"

