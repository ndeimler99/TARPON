#!/usr/bin/env python3

import gzip
import argparse
import pysam
import pandas as pd

def main(args):

    header = ""
    # stats_dict = {}

    stats_dict = pd.read_table(args.stats_fh, delimiter="\t")

    # with open(args.stats_fh, "r") as stats:
    #     linecount = 1
    #     for line in stats:
    #         if linecount == 1:
    #             linecount += 1
    #             header = line.strip()
    #             continue
    #         line = line.strip().split()
    #         stats_dict[line[0]] = line
    
    # identify unique alignments from aln file
    aln_file = pysam.AlignmentFile(args.alignment, "r", check_sq=False)
    alignment_dict = {}
    for aln in aln_file:
        if aln.query_name in alignment_dict:
            alignment_dict[aln.query_name].append(aln.reference_name)
        else:
            alignment_dict[aln.query_name] = [aln.reference_name]

    with open(args.new_stats_fh, "w") as stats_fh:
        stats_fh.write("read_id\tCluster\n")
        for read in stats_dict["read_id"]:
            if read in alignment_dict:
                if len(alignment_dict[read]) == 1:
                    stats_fh.write("{}\t{}\n".format(read, alignment_dict[read][0]))
                else:
                    stats_fh.write("{}\tNA\n".format(read))
            
            else:
                stats_fh.write("{}\tNA\n".format(read))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_fh", required=True)
    parser.add_argument("--new_stats_fh", required=True)
    parser.add_argument("--alignment", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
