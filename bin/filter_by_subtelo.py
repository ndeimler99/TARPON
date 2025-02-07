#!/usr/bin/env python3

import argparse
import regex
import pysam


def main(args):

    # convert parameters from strings to usable data types
    args.min_subtelo_length = int(args.min_subtelo_length)
    args.min_subtelo_threshold = float(args.min_subtelo_threshold)
    
    fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    subtelo_pass = pysam.AlignmentFile(args.passes_subtelo, "wb", template=fh)
    subtelo_fail = pysam.AlignmentFile(args.fails_subtelo, "wb", template=fh)

    if args.mutant == "false":
        pass
    else:
        hamming = 0
        for i in range(0, len(args.mutant)):
            if args.mutant[i] != args.repeat[i]:
                hamming += 1
        if hamming <= 1:
            args.mutant = "false"

    for aln in fh:    
        # ensure that the read length is not less than the length of minimum subtelomere stretch
        if len(aln.query_sequence) < args.min_subtelo_threshold:
            subtelo_fail.write(aln)
        else:
            matches = list(regex.finditer(r'(%s){s<=1}' % args.repeat, aln.query_sequence[0:args.min_subtelo_length]))
                    # if the read is long enough check the frequency of one nucleotide deviation of the telomeric repeat
                    # if that frequency is less than the min_subtelo_threshold the read passes
                    # if it is greater than min_subtelo_threshold the read begins in subtelomeric DNA and is excluded from further analysis
            if args.mutant == "false":
                if len(matches) * len(args.repeat) / args.min_subtelo_length <= args.min_subtelo_threshold:
                    subtelo_pass.write(aln)
                else:
                    subtelo_fail.write(aln)
            else:
                if (len(matches) * len(args.repeat) + aln.query_sequence[0:args.min_subtelo_length].count(args.mutant) * len(args.mutant)) / args.min_subtelo_length <= args.min_subtelo_threshold:
                    subtelo_pass.write(aln)
                else:
                    subtelo_fail.write(aln)

    fh.close()
    subtelo_pass.close()
    subtelo_fail.close()



def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--min_subtelo_length", required=True)
    parser.add_argument("--min_subtelo_threshold", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--passes_subtelo", required=True)
    parser.add_argument("--fails_subtelo", required=True)
    parser.add_argument("--mutant", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)