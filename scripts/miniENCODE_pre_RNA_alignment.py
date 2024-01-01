#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes RNA-seq data using HISAT2.

import sys
import getopt
import re
import subprocess
import os
import logging
import time
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

def rna_seq_gz(dataDir, hisat2, ht2index, bed, n, sra, gz, trim):
    fq_dir = f'{dataDir}/fastq/'
    align_dir = f'{dataDir}/alignment/'
    layout = sra[1]
    sra = sra[0]
    sra_bam = f'{align_dir}{sra}.bam'
    n = str(n)
    gzsuffix = '.gz' if gz == "True" else ''

    if layout == 'paired':
        fq1 = f'{fq_dir}{sra}_1_val_1.fq{gzsuffix}'
        fq2 = f'{fq_dir}{sra}_2_val_2.fq{gzsuffix}'
        mf.run_shell(f'{hisat2} -p {n} -x {ht2index} -1 {fq1} -2 {fq2} | samtools sort -@ {n} -o {sra_bam}')
    elif layout == 'single':
        fq = f'{fq_dir}{sra}_trimmed.fq{gzsuffix}'
        mf.run_shell(f'{hisat2} -p {n} -x {ht2index} -U {fq} | samtools sort -@ {n} -o {sra_bam}')

ht2index = f'{outpathway}/reference/index_hisat2/GRCz11'
bed = f'{outpathway}/reference/GRCz11_chr1_100.bed'
hisat2 = 'hisat2'

mf.run_shell(f'mkdir -p {outpathway}/alignment')

with open(inputfile, 'r') as SraListFile:
    for line in SraListFile:
        line = line.strip()
        log.info(line)
        line = line.split(',')
        rna_seq_gz(outpathway, hisat2, ht2index, bed, threads, line, gz='True', trim='True')
