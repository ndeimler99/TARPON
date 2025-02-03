#!/usr/bin/env python3

import argparse
import regex
import matplotlib.pyplot as plt
import pysam


def main(args):

    # convert parameters from strings to usable data types

    args.sliding_window = int(args.sliding_window)
    args.sliding_window_interval = int(args.sliding_window_interval)

    # create telomere statistics dictionary to draw telo start
    start_dict = {}
    linecount = 0
    with open(args.telo_stats, 'r') as telo_stats_fh:
        for line in telo_stats_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            start_dict[line[0]] = int(line[3])

    # for each read in fastq file draw the telomeric sequence
    input_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    for aln in input_fh:
        fig1, ax1 = plt.subplots(1,1)
        fig1.set_size_inches(12,4)
        ax1.set_xlabel("Sliding Window Start Position")
        ax1.set_ylabel("Percentage of Sliding\nWindow that is Telomeric")
        ax1.set_ylim(0,100)
        ax1.vlines(start_dict[aln.query_name], 0, 100, label="VRR Start", color="red")
        ax1.set_title('{}-{} Strand'.format(aln.query_name, aln.get_tag("XS")))
        x = []
        perf = []
        subs = []
        for i in range(0, len(aln.query_sequence)-args.sliding_window, args.sliding_window_interval):
            x.append(i)
            perf.append((aln.query_sequence[i:i+args.sliding_window].count(args.repeat) * len(args.repeat) / args.sliding_window) * 100)
            subs.append((len(list(regex.finditer(r"(%s){s<=1}" % args.repeat, aln.query_sequence[i:i+args.sliding_window]))) * len(args.repeat)) / args.sliding_window * 100)
        
        ax1.plot(x, perf, label="Perfect Repeats")
        ax1.plot(x, subs, label="One Nucl. Substitution")
        ax1.legend()
        #print(read[0].split()[0].strip("@"))
        fig1.savefig("{}.pdf".format(aln.query_name.strip("@")), format="pdf")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--telo_stats", required=True)
    parser.add_argument("--sliding_window", required=True)
    parser.add_argument("--sliding_window_interval", required=True)

    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)