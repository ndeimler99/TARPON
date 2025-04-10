#!/usr/bin/env python3

import gzip
import argparse
import pysam

def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def main(args):

    # convert parameters from strings to usable data types

    args.c_strand_only = args.c_strand_only == "true"
    args.repeat_count = int(args.repeat_count)

    # for read in gzipped fastq file perform analysis
    # if the argument c strand only is set convert the telomeric repeat is C strand and identify reads
    # otherwise check for both forward and reverse telomeric repeats
    input_file_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    out_fh = pysam.AlignmentFile(args.out_file, "wb", template=input_file_fh)
    non_telo_fh = pysam.AlignmentFile(args.non_telo, "wb", template=input_file_fh)

    for aln in input_file_fh:
        aln.query_name = aln.query_name.split()[0]
        if args.mutant == "false":
            if args.c_strand_only and aln.query_sequence.count(rev_complement(args.repeat)) >= args.repeat_count:
                out_fh.write(aln)
            elif aln.query_sequence.count(args.repeat) >= args.repeat_count or aln.query_sequence.count(rev_complement(args.repeat)) >= args.repeat_count:
                out_fh.write(aln)
            else:
                non_telo_fh.write(aln)
        else:
            if args.c_strand_only and aln.query_sequence.count(rev_complement(args.repeat)) + aln.query_sequence.count(rev_complement(args.mutant)) >= args.repeat_count:
                out_fh.write(aln)
            elif aln.query_sequence.count(args.repeat) + aln.query_sequence.count(args.mutant) >= args.repeat_count or aln.query_sequence.count(rev_complement(args.repeat)) + aln.query_sequence.count(rev_complement(args.mutant)) >= args.repeat_count:
                out_fh.write(aln)
            else:
                non_telo_fh.write(aln)
        
    input_file_fh.close()
    out_fh.close()
    non_telo_fh.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",required=True)
    parser.add_argument("--repeat",required=True)
    parser.add_argument("--repeat_count", required=True)
    parser.add_argument("--c_strand_only", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--non_telo", required=True)
    parser.add_argument("--mutant", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)