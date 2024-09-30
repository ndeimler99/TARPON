#!/usr/bin/env python3

import argparse
import regex
import matplotlib.pyplot as plt


def main(args):
    start_dict = {}
    linecount = 0
    with open(telo_stats, 'r') as telo_stats_fh:
        for line in telo_stats_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            start_dict[line[0]] = int(line[3])

    with open(input_file, "r") as input_fh:
        linecount = 0
        read = []
        for line in input_fh:

            linecount += 1
            read.append(line.strip())
            if linecount % 4 == 0:
                print(read[0])
                fig1, ax1 = plt.subplots(1,1)
                fig1.set_size_inches(12,4)
                ax1.set_xlabel("Sliding Window Start Position")
                ax1.set_ylabel("Percentage of Sliding\nWindow that is Telomeric")
                ax1.set_ylim(0,100)
                ax1.vlines(start_dict[read[0].split()[0]], 0, 100, label="VRR Start", color="red")
                ax1.set_title('{}-{} Strand'.format(read[0].split()[0], read[0].split()[1]))
                x = []
                perf = []
                subs = []
                for i in range(0, len(read[1])-sliding_window, sliding_window_interval):
                    x.append(i)
                    perf.append((read[1][i:i+sliding_window].count(repeat) * len(repeat) / sliding_window) * 100)
                    subs.append((len(list(regex.finditer(r"(%s){s<=1}" % repeat, read[1][i:i+sliding_window]))) * len(repeat)) / sliding_window * 100)
                
                ax1.plot(x, perf, label="Perfect Repeats")
                ax1.plot(x, subs, label="One Nucl. Substitution")
                ax1.legend()
                #print(read[0].split()[0].strip("@"))
                fig1.savefig("{}.pdf".format(read[0].split()[0].strip("@")), format="pdf")
                read = []

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