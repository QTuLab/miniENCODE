#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: Read a matrix from a file and output the maximum Z-score for each row.

import sys

# Open the matrix file provided as a command line argument
with open(sys.argv[1], 'r') as matrix:
    for line in matrix:
        line = line.rstrip().split("\t")
        rDHS = line[0]
        ZscoreArray = [float(i) for i in line[1:]]
        print(f"{rDHS}\t{max(ZscoreArray)}")
