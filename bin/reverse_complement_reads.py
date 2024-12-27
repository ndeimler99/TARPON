#!/usr/bin/env python3

import argparse

def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def main(args):
    args.c_strand_only = args.c_strand_only == "true"
    args.threshold = float(args.threshold)

    with open(args.input_file, 'r') as input_file_fh, open(args.out_file, 'w') as out_fh, open(args.removed_reads, 'w') as filtered:
        linecount = 0
        read = []
        for line in input_file_fh:
            linecount += 1
            read.append(line.strip())
            if linecount % 4 == 0:
                c = read[1].count(rev_complement(args.repeat))
                g = read[1].count(args.repeat)
                if args.c_strand_only:
                    if c/(c+g) >= args.threshold:
                        # write out to main file and c_strand file
                        out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                    else:
                        filtered.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'None', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                else:
                    if c/(c+g) <= args.threshold and c/(c+g) >= 1-args.threshold:
                        filtered.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'None', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                    elif c/(c+g) >= args.threshold:
                        #write out to main file and c_strand file
                        out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                    elif g/(c+g) >= args.threshold:
                        #write out to main file and g_strand file
                        out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'G', read[1], read[2], read[3]))
                    
                read = []

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--threshold", required=True)
    parser.add_argument("--c_strand_only", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--removed_reads", required=True)
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)