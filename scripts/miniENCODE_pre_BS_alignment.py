#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script aligns BS-seq data using Bismark.

import sys
import getopt 
import time
import subprocess
import chardet
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

def b2s(byte):
    encode_type = chardet.detect(byte)
    byte = byte.decode(encode_type['encoding'])
    return byte

# Define function to check layout of SRA data
def layout_check(sra):
    srafile = 'SRA/' + line + '.sra'
    lay = mf.run_shell('fastq-dump -X 1 --split-spot -Z ' + srafile + ' | wc -l')
    if (b2s(lay.stdout).strip()=='8'):
        layout = 'paired'
    elif b2s(lay.stdout).strip()=='4':
        layout = 'single'
    return layout

# Define function for methylation data
def bs_seq_gz(dataDir, bisindex, n, sra):
    fq_dir = f'{dataDir}/fastq/'
    align_dir = f'{dataDir}/alignment/'
    layout = sra[1]
    sra = sra[0]
    n = str(n)
    if layout == 'paired':
        fq1 = f'{fq_dir}{sra}_1_val_1.fq.gz'
        fq2 = f'{fq_dir}{sra}_2_val_2.fq.gz'
        mf.run_shell(f'bismark --parallel {n} -o {align_dir} --genome {bisindex} -1 {fq1} -2 {fq2}')
        mf.run_shell(f'deduplicate_bismark -o {sra} --output_dir {align_dir} --bam {align_dir}{sra}_1_val_1_bismark_bt2_pe.bam')
    elif layout == 'single':
        mf.run_shell(f'bismark --parallel {n} -o {align_dir} --genome {bisindex} {fq_dir}{sra}_trimmed.fq.gz')
        mf.run_shell(f'deduplicate_bismark -o {sra} --output_dir {align_dir} --bam {align_dir}{sra}_trimmed_bismark_bt2.bam')

mf.run_shell(f'mkdir -p {outpathway}/alignment')

# Bismark index directory
bisindex = f'{outpathway}/reference/index_bismark/'

# Process from the input file
SraListFile = open(inputfile, 'r')
for line in SraListFile:
    line = line.strip()
    log.info(line)
    line = line.split(',')
    bs_seq_gz(outpathway, bisindex, 15, line)
    print()
SraListFile.close()
