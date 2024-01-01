#!/bin/bash

# miniENCODE preprocessing pipeline
# Script Purpose: Process ELS (enhancer-like signature) data for ATAC-seq and H3K27ac ChIP-seq.

# Define the data directory
dir=$PWD
dataDir="${dir}/data/ELS/"
peaks="${dataDir}/rDAR_final.bed"
scriptDir="${dir}/scripts"
output="${dataDir}/process"

# Set Python environment
pyenv global system

# Loop through sample and mode combinations
for sample in Blastula Blood
do
    for mode in ATAC H3K27ac
    do
        echo "Processing $sample in $mode mode"
        width=0
        bw="$dataDir/$sample/${mode}_${sample}.bam.bw"
        [ "$mode" == "H3K27ac" ] && width=500

        mkdir -p "${output}/${sample}/${mode}/"
        cd "${output}/${sample}/${mode}/" || exit

        # Generate a modified BED file ('little') with adjusted coordinates
        awk -F "\t" -v width="$width" '{printf "%s\t%.0f\t%.0f\t%s\n", $1, $2-width, $3+width, $4}' "$peaks" | \
        awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4; else print $0}' | \
        sort -u > little

        bigWigAverageOverBed -bedOut=out2.bed "$bw" little out2

        python "${scriptDir}/miniENCODE_ELS_log_zscore.py" out2 > l

        sort -k2,2rg l | \
        awk 'BEGIN {rank=0; before=0; running=1}{if ($2 != before) rank = running; \
        print $1 "\t" $2 "\t" $3 "\t" rank; before=$2; running += 1}' | \
        sort -k1,1 > out.txt

        # Create a matrix from the results
        echo "Creating Matrix..."
        paste out.txt | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=4) printf "%s\t", $i; print ""}' > matrix
        
        echo "Determining maxZ..."
        python "${scriptDir}/miniENCODE_ELS_max_zscore.py" matrix > maxZ.txt
    done
done

