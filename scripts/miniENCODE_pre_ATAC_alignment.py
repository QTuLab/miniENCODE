#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script processes ATAC-seq data using Bowtie2 for alignment
# and Picard for marking duplicates.

import sys
import getopt
import chardet
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
    log.error("Usage: " + sys.argv[0] + " -i|--input <inputSRAListfile> \
              -o|--output <outpathway> -t|--threads <numberOfThreads>")
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:', ['input=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

inputfile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-i', '--input'):
        inputfile = arg
    if opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if ((inputfile == '') or (outpathway == '')):
    usage()
    sys.exit(2)

# Define function to convert bytes to string using chardet
def b2s(byte):
    encode_type = chardet.detect(byte)
    byte = byte.decode(encode_type['encoding'])
    return byte

# Define function for ATAC-seq processing
def atac_bt2(dataDir, bt2base, n, sra):
    fq_dir = f'{dataDir}/fastq/'
    align_dir = f'{dataDir}/alignment/'
    layout = sra[1]
    sra = sra[0]
    n = str(n)
    mapbam = f'{align_dir}{sra}.map.bam'
    dupstxt = f'{align_dir}{sra}.dups.txt'
    deduplen = f'{align_dir}{sra}.dedup.length.txt'
    finalbam = f'{align_dir}{sra}.final.bam'
    log.info('Running Bowtie2')
    if layout == 'paired':
        fq1 = f'{fq_dir}{sra}_1_val_1.fq.gz'
        fq2 = f'{fq_dir}{sra}_2_val_2.fq.gz'
        mf.run_shell(f'bowtie2 -X 2000 --qc-filter --mm -p {n} -x {bt2base} -1 {fq1} -2 {fq2} | samtools view -@ {n} -1 -S - | samtools sort -@ {n} - -o {mapbam}')
    elif layout == 'single':
        fq = f'{fq_dir}{sra}_trimmed.fq.gz'
        mf.run_shell(f'bowtie2 -X 2000 --qc-filter --mm -p {n} -x {bt2base} -U {fq} | samtools view -@ {n} -1 -S - | samtools sort -@ {n} - -o {mapbam}')
    log.info('Running Picard')
    mf.run_shell(f'picard-tools MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true I={mapbam} O={finalbam} M={dupstxt}')
    log.info('Calculating deduplicated read lengths')
    mf.run_shell(f'samtools view -O SAM -q 30 -F 4 -@ {n} {finalbam} | cut -f9 | sed "s/-//g" | sort | uniq -c > {deduplen}')

############################
mf.run_shell(f'mkdir -p {outpathway}/alignment')

# Bowtie2 index directory
bt2base = f'{outpathway}/reference/index_bowtie2/GRCz11'

# Process from the input file
with open(inputfile, 'r') as SraListFile:
    for line in SraListFile:
        line = line.strip()
        log.info(line)
        line = line.split(',')
        atac_bt2(outpathway, bt2base, threads, line)
