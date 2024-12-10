#!/usr/bin/env python3

import argparse
import regex


def main(args):

    # convert parameters from strings to usable data types
    args.min_subtelo_length = int(args.min_subtelo_length)
    args.min_subtelo_threshold = float(args.min_subtelo_threshold)
    
    with open(args.input_file, 'r') as fh, \
    open(args.passes_subtelo, 'w') as subtelo_pass, open(args.fails_subtelo, 'w') as subtelo_fail:
        read = []
        linecount = 0
        # for every telomeric sequence in the fastq file
        for line in fh:
            linecount += 1
            if linecount % 4 == 0:
                read.append(line.strip())
                # ensure that the read length is not less than the length of minimum subtelomere stretch
                if len(read[1]) < args.min_subtelo_threshold:
                    subtelo_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                else:
                    matches = list(regex.finditer(r'(%s){s<=1}' % args.repeat, read[1][0:args.min_subtelo_length]))
                    # if the read is long enough check the frequency of one nucleotide deviation of the telomeric repeat
                    # if that frequency is less than the min_subtelo_threshold the read passes
                    # if it is greater than min_subtelo_threshold the read begins in subtelomeric DNA and is excluded from further analysis
                    if len(matches) * len(args.repeat)/ args.min_subtelo_length <= args.min_subtelo_threshold:
                        subtelo_pass.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                    else:
                        subtelo_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                read = []
            else:
                read.append(line.strip())

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--min_subtelo_length", required=True)
    parser.add_argument("--min_subtelo_threshold", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--passes_subtelo", required=True)
    parser.add_argument("--fails_subtelo", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)