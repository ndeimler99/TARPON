#!/usr/bin/env python3

import gzip
import argparse


def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def main(args):
    args.c_strand_only = args.c_strand_only == "true"
    args.repeat_count = int(args.repeat_count)

    with gzip.open(args.input_file, 'rt') as input_file_fh, open(args.out_file, 'w') as out_fh, open(args.non_telo, "w") as non_telo_fh:
        linecount = 0
        read = []
        for line in input_file_fh:
            linecount += 1
            read.append(line.strip())
            if linecount % 4 == 0:
                if args.c_strand_only and read[1].count(rev_complement(args.repeat)) >= args.repeat_count:
                    out_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3])) 
                elif read[1].count(args.repeat) >= args.repeat_count or read[1].count(rev_complement(args.repeat)) >= args.repeat_count:
                    out_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3]))    
                else:
                    non_telo_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3]))
                read = []

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",required=True)
    parser.add_argument("--repeat",required=True)
    parser.add_argument("--repeat_count", required=True)
    parser.add_argument("--c_strand_only", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--non_telo", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)