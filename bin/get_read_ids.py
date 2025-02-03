#!/usr/bin/env python3

import gzip
import argparse
import pysam

def main(args):
    input_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    with open(args.output_file, "w") as out_fh:
        for aln in input_fh:
            out_fh.write(aln.query_name + "\n")
    input_fh.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)