#!/usr/bin/env python3

import argparse
import pysam

def main(args):
    input_file_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    g_fh = pysam.AlignmentFile(args.g_file, "wb", template=input_file_fh)
    c_fh = pysam.AlignmentFile(args.c_file, "wb", template=input_file_fh)

    for aln in input_file_fh:
        if aln.get_tag("XS") == "G":
            g_fh.write(aln)
        else:
            c_fh.write(aln)
    
    input_file_fh.close()
    g_fh.close()
    c_fh.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--g_file", required=True)
    parser.add_argument("--c_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)