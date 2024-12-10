#!/usr/bin/env python3

import argparse
import sys

def main(args):
    # Opens sample file and creates a list of barcodes to check
    barcodes = []
    with open(args.sample_file, 'r') as sample_fh:
        linecount = 0
        for line in sample_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split(",")
            barcodes.append(line[1])

    # iterates through list of barcodes and compares each barcode to all other barcodes
    # checking that a barcode does not have a minimum hamming distance less than the 
    # number of allowed errors used during it's identification in a telomeric sequence
    for i in range(0, len(barcodes)-1):
        for j in range(i+1, len(barcodes)):
            if barcodes[i] == barcodes[j]:
                sys.exit(1)
            else:
                distance = 0
                for k in range(0, len(barcodes[i])):
                    if barcodes[i][k] != barcodes[j][k]:
                        distance += 1
                if distance <= int(args.barcode_errors):
                    sys.exit(1)

def argparser():

    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_file", required=True, help="comma separated Sample file with first column as sample Id and second column as barcode")
    parser.add_argument("--barcode_errors", required=True, help="Number of errors allowed in a barcode to be identified in a telomeric sequence")
    return parser
    
if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)


        
        
