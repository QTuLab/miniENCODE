#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: This script merges ChIP-seq data and calls peaks using MACS2.

import sys
import getopt
import time
import subprocess
import os
import logging
import signal

# Configure logging to show info level messages
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Read command line arguments
def usage():
    print("Error!")
    print("Usage: " + sys.argv[0] + " -s|--sample <sampleListFile>  \
            -c|--callpeak <callpeakFile> -o|--output <outpathway> \
            -t|--threads <numberOfThreads>")
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 's:c:o:t:', [
                               'sample=', 'callpeak=', 'output=', 'threads='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

samplefile = ''
callpeakfile = ''
outpathway = ''
threads = 1

for opt, arg in opts:
    if opt in ('-s', '--sample'):
        samplefile = arg
    if opt in ('-c', '--callpeak'):
        callpeakfile = arg
    if opt in ('-o', '--output'):
        outpathway = arg
    elif opt in ('-t', '--threads'):
        threads = arg

if ((samplefile == '') or (callpeakfile == '') or (outpathway == '')):
    usage()
    sys.exit(2)

def get_ticks():
    return getattr(time, 'perf_counter', getattr(time, 'time'))()

# Define function to run shell commands with error handling
def run_shell(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    t0 = get_ticks()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    t1 = get_ticks()
    err_str = ('PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}, STDERR={stde}, STDOUT={stdo}').format(
        pid=pid, pgid=pgid, rc=rc, dur=t1 - t0, stde=stderr.strip(), stdo=stdout.strip())
    if rc:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')

# Function to merge ChIP-seq data
def chip_seq_merge(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    bam_file = f'{align_dir}{relist[0]}.final.bam'
    n = str(n)
    if len(relist) == 2:
        log.info('Rename data')
        run_shell(f'mv {align_dir}{relist[1]}.final.bam {bam_file}')
    elif len(relist) > 2:
        nfile = len(relist)
        log.info(f'Merge {nfile} files')
        merge_cmd = ' '.join([f'{align_dir}{file}.final.bam' for file in relist])
        run_shell(f'samtools merge -f -@ {n} {merge_cmd}')
    run_shell(f'samtools sort -@ {n} -o {bam_file} {bam_file}')
    run_shell(f'samtools index -@ {n} {bam_file}')
    run_shell(f'bamCoverage -bs 1 --normalizeUsing RPKM --numberOfProcessors {n} -b {bam_file} -o {bam_file}.bw')

# Peak calling
def chip_callpeak(dataDir, relist, n):
    align_dir = f'{dataDir}/alignment/'
    peak_dir = f'{dataDir}/peak/'
    n = str(n)
    log.info('Peak calling')
    if len(relist) == 1:
        run_shell(f'macs2 callpeak --outdir {peak_dir} --keep-dup all -g 1.1e9 -t {align_dir}{relist[0]}.final.bam -n {relist[0]}')
    elif len(relist) == 2:
        run_shell(f'macs2 callpeak --outdir {peak_dir} --keep-dup all -g 1.1e9 -t {align_dir}{relist[0]}.final.bam -c {align_dir}{relist[1]}.final.bam -n {relist[0]}')

###################################

# Process data from the sample file
with open(samplefile, 'r') as SampleListFile:
    for line in SampleListFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        chip_seq_merge(outpathway, relist, threads)

# Process data to call peaks
with open(callpeakfile, 'r') as CallPeakFile:
    for line in CallPeakFile:
        line = line.strip()
        log.info(line)
        relist = line.split(',')
        chip_callpeak(outpathway, relist, threads)
