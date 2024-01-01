#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script quantifies RNA-seq data using StringTie.

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

def rna_quantify(dataDir, relist, n, gff3):
    align_dir = f'{dataDir}/alignment/'
    name = relist[0]
    bam_name = f'{align_dir}{name}.bam'
    qt_dir = f'{dataDir}/quant/'
    mf.run_shell(f'stringtie -G {gff3} -p {n} -e -B -A {qt_dir}{name}_gene_abundance.tab -C {qt_dir}{name}_cov_refs.gtf -l {name} -o {qt_dir}{name}.gtf {bam_name}')

gff3 = f'{outpathway}/reference/GRCz11_chr1_100.gff3'
mf.run_shell(f'mkdir -p {outpathway}/quant')

with open(samplefile, 'r') as SampleListFile:
    for line in SampleListFile:
        line = line.strip()
        relist = line.split(',')
        rna_quantify(outpathway, relist, threads, gff3)
