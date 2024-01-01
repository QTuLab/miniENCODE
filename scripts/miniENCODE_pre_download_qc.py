#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script performs quality control on SRA data and converts it to fastq format.

import sys
import getopt
import subprocess
import time
import re
import chardet
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

def b2s(byte):
    return byte.decode(chardet.detect(byte)['encoding'])

def check_qc(dataDir, sra, n):
    sra_dir = f'{dataDir}/srafile/'
    fq_dir = f'{dataDir}/fastq/'
    qc_dir = f'{dataDir}/qc/'
    layout = sra[1]
    sra = sra[0]
    srafile = f'{sra_dir}/{sra}/{sra}.sra'
    log.info('Downloading SRA file')
    mf.run_shell(f'prefetch -p -O {sra_dir} {sra}')
    vdb = subprocess.run(f'vdb-validate {srafile}', shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    vdbcheck = 'is consistent' in b2s(vdb.stderr)
    if vdbcheck:
        log.info('SRA check successful')
        log.info('Starting SRA to fastq conversion and quality control')
        if layout == 'paired':
            log.info('Layout: paired')
            mf.run_shell(f'fastq-dump --gzip -O {fq_dir} --split-files {srafile}')
            mf.run_shell(f'fastqc -t {n} -o {qc_dir} {fq_dir}/{sra}_*.fastq.gz')
        elif layout == 'single':
            log.info('Layout: single')
            mf.run_shell(f'fastq-dump --gzip -O {fq_dir} {srafile}')
            mf.run_shell(f'fastqc -t {n} -o {qc_dir} {fq_dir}/{sra}.fastq.gz')
    else:
        log.error('Error in SRA check')

# Create necessary directories
mf.run_shell(f'mkdir -p {outpathway}/srafile {outpathway}/fastq {outpathway}/qc')

# Read input file line by line and process data
with open(inputfile, 'r') as SraListFile:
    for line in SraListFile:
        line = line.strip()
        log.info(line)
        sra_info = line.split(',')
        check_qc(outpathway, sra_info, threads)
