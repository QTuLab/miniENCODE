#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: Process signal data, compute logarithmic values, and standardize based on mean and standard deviation.

import numpy as np
import sys
import math

# Check if input file argument is provided
if len(sys.argv) < 2:
    print("Usage: script.py <input_file>")
    sys.exit(1)

# Initialize lists to store signal data, master peak names, and calculations
sig = [[], []]
masterPeak = []
calculate = []

# Process the input file
with open(sys.argv[1], 'r') as bigWig:
    for line in bigWig:
        line = line.rstrip().split("\t")
        value = float(line[4])

        # Handle zero and non-zero values differently
        if value == 0:
            sig[1].append("Zero")
            sig[0].append(value)
            masterPeak.append(line[0])
        else:
            logged_value = math.log(value, 10)
            sig[1].append(logged_value)
            sig[0].append(value)
            calculate.append(logged_value)
            masterPeak.append(line[0])

# Calculate the mean and standard deviation
lmean = np.mean(calculate)
lstd = np.std(calculate)

# Output the processed data
for i, entry in enumerate(sig[1]):
    if entry != "Zero":
        normalized_value = (entry - lmean) / lstd
        print(f"{masterPeak[i]}\t{normalized_value}\t{sig[0][i]}\t{sig[1][i]}")
    else:
        print(f"{masterPeak[i]}\t{-10}\t{0}\t{0}")
