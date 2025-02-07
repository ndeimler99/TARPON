#!/usr/bin/env python3

import argparse
import pysam

def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def main(args):
    args.c_strand_only = args.c_strand_only == "true"
    args.threshold = float(args.threshold)

    input_file_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    out_fh = pysam.AlignmentFile(args.out_file, "wb", template=input_file_fh)
    filtered = pysam.AlignmentFile(args.removed_reads, "wb", template=input_file_fh)

    for aln in input_file_fh:
        if args.mutant == "false":
            c = aln.query_sequence.count(rev_complement(args.repeat))
            g = aln.query_sequence.count(args.repeat)
        else:
            c = aln.query_sequence.count(rev_complement(args.repeat)) + aln.query_sequence.count(rev_complement(args.mutant))
            g = aln.query_sequence.count(args.repeat) + aln.query_sequence.count(args.mutant)
        if args.c_strand_only:
            if c/(c+g) >= args.threshold:
                q = aln.query_qualities
                aln.query_sequence = rev_complement(aln.query_sequence)
                aln.set_tag("XS", "C")
                aln.query_qualities = q[::-1]
                out_fh.write(aln)
            else:
                filtered.write(aln)
        else:
            if c/(c+g) <= args.threshold and c/(c+g) >= 1-args.threshold:
                filtered.write(aln)
            elif c/(c+g) >= args.threshold:
                aln.set_tag("XS", "C")
                q = aln.query_qualities
                aln.query_sequence = rev_complement(aln.query_sequence)
                aln.query_qualities = q[::-1]
                out_fh.write(aln)
            elif g/(c+g) >= args.threshold:
                aln.set_tag("XS", "G")
                out_fh.write(aln)

    input_file_fh.close()
    out_fh.close()
    filtered.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--threshold", required=True)
    parser.add_argument("--c_strand_only", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--removed_reads", required=True)
    parser.add_argument("--mutant", required=True)
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)