#!/usr/bin/env python3

import gzip
import argparse

def main(args):
    with open(args.input_file, "r") as input_fh, open(args.output_file, "w") as out_fh:
        linecount  = 0 
        for line in input_fh: 
            linecount += 1
            if linecount % 4 == 1:
                line = line.strip().split()[0].strip("@")
                out_fh.write(line + "\n")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)