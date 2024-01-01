#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script calls peaks in ATAC-seq data using MACS2.

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
    log.error("Usage: " +
          sys.argv[0] + " -s|--sample <sampleListfile> \
            -o|--output <outpathway> -t|--threads <numberOfThreads>")
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 's:o:t:', [
                               'sample=', 'output=', 'threads='])
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

def atac_seq_merge(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    peak_dir = f'{dataDir}/peak/'
    bam_file = f'{align_dir}{relist[0]}.final.bam'
    n = str(n)
    if len(relist) == 2:
        log.info('Rename data')
        mf.run_shell(f'mv {align_dir}{relist[1]}.final.bam {bam_file}')
    elif len(relist) > 2:
        nfile = len(relist)
        log.info(f'Merge {nfile} files')
        merge_cmd = ' '.join([f'{align_dir}{file}.final.bam' for file in relist])
        mf.run_shell(f'samtools merge -f -@ {n} {merge_cmd}')
    mf.run_shell(f'samtools sort -@ {n} -o {bam_file} {bam_file}')
    mf.run_shell(f'samtools index -@ {n} {bam_file}')
    mf.run_shell(f'bamCoverage -bs 1 --normalizeUsing RPKM --numberOfProcessors {n} -b {bam_file} -o {bam_file}.bw')
    log.info('Peak calling')
    mf.run_shell(f'macs2 callpeak --outdir {peak_dir} --keep-dup all -g 1.1e9 -t {bam_file} -n {relist[0]}')

# Process data from the sample file
with open(samplefile, 'r') as SampleListFile:
    for line in SampleListFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        atac_seq_merge(outpathway, relist, threads)
