#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script merges RNA-seq data and generates BigWig files.

import sys
import getopt
import time
import subprocess
import os
import logging
import signal
import miniENCODE_function as mf

logging.basicConfig(format='[%(asctime)s %(levelname)s] %(message)s', stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def usage():
    log.error("Error!")
    log.error("Usage: " + sys.argv[0] + " -s|--sample <sampleListfile> \
            -o|--output <outpathway> -t|--threads <numberOfThreads>")
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], 's:o:t:', ['sample=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()

samplefile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-s', '--sample'):
        samplefile = arg
    elif opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if not samplefile or not outpathway:
    usage()

def rna_seq_merge(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    name = relist[0]
    bam_name = f'{align_dir}{name}.bam'
    n = str(n)

    if len(relist) == 2:
        mf.run_shell(f'mv {align_dir}{relist[1]}.bam {bam_name}')
    elif len(relist) > 2:
        log.info('Merging data')
        merge_cmd = ' '.join([f'{align_dir}{file}.bam' for file in relist])
        mf.run_shell(f'samtools merge -f -@ {n} {merge_cmd}')
    mf.run_shell(f'samtools sort -@ {n} -o {bam_name} {bam_name}')
    mf.run_shell(f'samtools index -@ {n} {bam_name}')
    mf.run_shell(f'bamCoverage -bs 1 --normalizeUsing RPKM --numberOfProcessors {n} -b {bam_name} -o {bam_name}.bw')

with open(samplefile, 'r') as SampleListFile:
    for line in SampleListFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        rna_seq_merge(outpathway, relist, threads)
