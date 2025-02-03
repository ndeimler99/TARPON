#!/usr/bin/env python3

import argparse
import numpy as np
import pysam

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_files", nargs="+", required=True)
    return parser

def get_mean_qual(seq):
    """Returns the mean phred converted quality score for a given seq"""
    return sum([i for i in seq]) / len(seq)

def main(args):
    print("Sample\tnum_seqs\tMin\tQ1\tQ2\tQ3\tMax\tMean")
    for file in args.fastq_files:
        num_seqs = 0
        quality_list = []
        linecount = 0

        fh = pysam.AlignmentFile(file, "rb", check_sq=False)

        for aln in fh:
            num_seqs += 1
            quality_list.append(get_mean_qual(aln.query_qualities))
        fh.close()

        if len(quality_list) > 0:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(file, num_seqs, min(quality_list), np.quantile(quality_list, 0.25), np.quantile(quality_list,0.50), np.quantile(quality_list, 0.75), max(quality_list), np.mean(quality_list)))
        else:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(file, num_seqs, 0, 0, 0, 0, 0, 0)) 

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

