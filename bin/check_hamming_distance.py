#!/usr/bin/env python3

import argparse

def main(args):
    barcodes = []
    with open(args.sample_file, 'r') as sample_fh:
        linecount = 0
        for line in sample_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split(",")
            barcodes.append(line[1])

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
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--barcode_errors", required=True)
    return parser
    
if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)


        
        
