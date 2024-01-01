#!/usr/bin/python3

# miniENCODE preprocessing pipeline
# Script Purpose: Calculate the distance between peaks and transcription start sites (TSS) in genomic data.

import sys

def calculate_diff(mode, line):
    """
    Calculate the difference between peak and TSS positions.

    Args:
    mode (str): Mode of operation, either 'specific' or 'agnostic'.
    line (list): A list representing a line from the input file, split into elements.

    Returns:
    int: The calculated distance.
    """
    peak = (int(line[1]) + int(line[2])) / 2
    tss_index = 7 if mode == "specific" else 5
    direction_index = 11 if mode == "specific" else 9
    tss = int(line[tss_index])
    diff = peak - tss if line[direction_index] == "-" else tss - peak
    return diff

def main():
    """
    Main function to read the input file and calculate distances based on the specified mode.
    """
    if len(sys.argv) < 3:
        print("Usage: script.py <input_file> <mode>")
        sys.exit(1)

    input_file = sys.argv[1]
    mode = sys.argv[2]

    if mode not in ["specific", "agnostic"]:
        print("Mode should be 'specific' or 'agnostic'.")
        sys.exit(1)

    ccreDict = {}
    try:
        with open(input_file) as f:
            for line in f:
                line = line.rstrip().split("\t")
                diff = calculate_diff(mode, line)
                gene_name = line[3]

                if gene_name not in ccreDict or abs(diff) < abs(ccreDict[gene_name]):
                    ccreDict[gene_name] = diff

        for entry in ccreDict:
            print(f"{entry}\t{ccreDict[entry]}")

    except FileNotFoundError:
        print(f"File not found: {input_file}")
        sys.exit(1)

if __name__ == "__main__":
    main()
