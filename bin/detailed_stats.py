#!/usr/bin/env python3

import regex
import argparse

def get_mean_qual(seq):
    """Returns the mean phred converted quality score for a given seq"""
    return sum([ord(i)-33 for i in seq]) / len(seq)

def main(args):
    with open(args.fastq_file, "r") as input_fh, open(args.stats_in, "r") as stats_in, open(args.stats_out, "w") as stats_out:
        linecount = 0
        stats_dict = {}
        # read in pre-existing telomeric statistics and store in a dictionary where the key is the read id and the value is the entire line
        for line in stats_in:  
            if linecount == 0:
                linecount += 1
                stats_out.write(line.strip() + "\tperc_variant\tperc_perfect\n")
            else:
                stats_dict[line.strip().split()[0]] = line.strip()
        
        # for every telomeric sequence in args.fastq_file add the detailed stats to the stored dictionary and write out to a new file
        linecount = 0
        read = []
        for line in input_fh:
            linecount += 1
            read.append(line.strip())
            if linecount % 4 == 0:
                telo_seq = read[1][int(stats_dict[read[0].split()[0]].split()[3]):]
                perfect_perc = telo_seq.count(args.repeat) * len(args.repeat) / len(telo_seq) * 100
                one_subs = len(list(regex.finditer(r'(%s){s<=1}' % args.repeat, telo_seq))) * len(args.repeat) / len(telo_seq) * 100
                stats_out.write('{}\t{}\t{}\t{}\t{}\n'.format(stats_dict[read[0].split()[0]], one_subs, perfect_perc))
                read = []
            



def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--stats_in", required=True)
    parser.add_argument("--stats_out", required=True) 
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

        
