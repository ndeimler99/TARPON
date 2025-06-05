#!/usr/bin/env python3

import gzip
import argparse
import pysam

def main(args):

    args.pre_vrr = int(args.pre_vrr)
    args.post_vrr = int(args.post_vrr)

    stats_dict = {}
    linecount = 1
    with open(args.stats_fh, "r") as stats:
        for line in stats:
            if linecount == 1:
                linecount += 1
                continue
            line = line.strip().split()
            stats_dict[line[0]] = int(line[3])
    

    telo_file = pysam.AlignmentFile(args.telomeric_bam, "rb", check_sq=False)

    with open(args.out_fasta, "w") as out_fh:
        for aln in telo_file:
            read_id = aln.query_name
            start = 0 if stats_dict[aln.query_name] - args.pre_vrr < 0 else stats_dict[aln.query_name] - args.pre_vrr
            end = len(aln.query_sequence) if  stats_dict[aln.query_name] + args.post_vrr > len(aln.query_sequence) else stats_dict[aln.query_name] + args.post_vrr

            sequence = aln.query_sequence[start:end]
            out_fh.write(">{}\n{}\n".format(read_id, sequence))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telomeric_bam", required=True)
    parser.add_argument("--stats_fh", required=True)
    parser.add_argument("--pre_vrr", required=True)
    parser.add_argument("--post_vrr", required=True)
    parser.add_argument("--out_fasta", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
