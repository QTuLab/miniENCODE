#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script trims RNA-seq data and performs quality control using FastQC.

import sys
import getopt
import subprocess
import time
import os
import logging
import signal
import miniENCODE_function as mf

logging.basicConfig(format='[%(asctime)s %(levelname)s] %(message)s', stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def usage():
    log.error("Error!")
    log.error("Usage: " + sys.argv[0] + " -i|--input <inputSRAListfile> -o|--output <outpathway> -t|--threads <numberOfThreads>")
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:', ['input=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()

inputfile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-i', '--input'):
        inputfile = arg
    elif opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if not inputfile or not outpathway:
    usage()

def trim_gz(dataDir, sra, n):
    fq_dir = f'{dataDir}/fastq/'
    qc_dir = f'{dataDir}/qc/'
    layout = sra[1]
    sra = sra[0]
    n = str(n)

    if layout == 'paired':
        fq_pe = ' '.join([f'{fq_dir}{sra}_{reads}.fastq.gz' for reads in [1, 2]])
        mf.run_shell(f'trim_galore --gzip --paired -j {n} -o {fq_dir} {fq_pe}')
    elif layout == 'single':
        fq_se = f'{fq_dir}{sra}.fastq.gz'
        mf.run_shell(f'trim_galore -j {n} -o {fq_dir} {fq_se}')

    mf.run_shell(f'fastqc -o {qc_dir} -t {n} {fq_dir}{sra}*.fq.gz')

# Create necessary directories
mk_dir = ' '.join([f'{outpathway}/{dirs}' for dirs in ["fastq", "qc"]])
mf.run_shell(f'mkdir -p {mk_dir}')

with open(inputfile, 'r') as SraListFile:
    for line in SraListFile:
        line = line.strip()
        line = line.split(',')
        trim_gz(outpathway, line, threads)
