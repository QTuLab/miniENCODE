#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes BS-seq data by renaming, merging, and converting files to bedGraph and BigWig formats.

import sys
import getopt
import time
import subprocess
import os
import logging
import signal
import miniENCODE_function as mf

# Configure logging to show info level messages
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Read command line arguments
def usage():
    log.error("Error!")
    log.error("Usage: " + sys.argv[0] + " -s|--sample <sampleListfile> \
              -o|--output <outpathway> -t|--threads <numberOfThreads>")
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 's:o:t:', ['sample=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

samplefile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-s', '--sample'):
        samplefile = arg
    if opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if ((samplefile == '') or (outpathway == '')):
    usage()
    sys.exit(2)

def bs_seq_bw(dataDir, ChromSizes, n, relist):
    align_dir = f'{dataDir}/alignment/'
    name_dir = f'{align_dir}/{relist[0]}'
    name = relist[0]
    n = str(n)

    if len(relist) == 2:
        log.info('Rename data')
        mf.run_shell(f'mv {align_dir}/{relist[1]}.deduplicated.bam {name_dir}.deduplicated.bam')

    elif len(relist) > 2:
        log.info('merge data')
        nfile = len(relist)
        log.info('Merge data {nfile} files')
        merge_cmd = ' '.join([f'{align_dir}{file}.deduplicated.bam' for file in relist])
        mf.run_shell(f'samtools merge -f -@ {n} {merge_cmd}')

    mf.run_shell(f'samtools sort -n -@ {n} -O BAM -o {name_dir}_sorted.bam {name_dir}.deduplicated.bam')
    log.info('bismark_methylation_extractor')
    mf.run_shell(f'bismark_methylation_extractor --gzip --bedGraph -o {align_dir} {name_dir}_sorted.bam')
    mf.run_shell(f'gzip -f -d {name_dir}_sorted.bedGraph.gz')
    subprocess.run(f"cat {name_dir}_sorted.bedGraph | sed '1d' | sort -k1,1 -k2,2n > {name_dir}_final.bedGraph", shell=True)
    mf.run_shell(f'bedClip {name_dir}_final.bedGraph {ChromSizes} {name_dir}_final_out.bedGraph')
    mf.run_shell(f'bedGraphToBigWig {name_dir}_final_out.bedGraph {ChromSizes} {name_dir}_final.bw')

# Chromosome sizes file
ChromSizes = f'{outpathway}/reference/GRCz11_chr1.chrom.sizes'

# Process data from the sample file
SraListFile = open(samplefile, 'r')
for line in SraListFile:
    line = line.strip()
    log.info(line)
    relist = line.split(',')
    bs_seq_bw(outpathway, ChromSizes, threads, relist)

SraListFile.close()
