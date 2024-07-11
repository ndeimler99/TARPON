#!/usr/bin/env python3

import sys
import regex

fastq_file = sys.argv[1]
repeat = sys.argv[2]
stats_in = sys.argv[3]
stats_out = sys.argv[4]

def get_mean_qual(seq):
    return sum([ord(i)-33 for i in seq]) / len(seq)

with open(fastq_file, "r") as input_fh, open(stats_in, "r") as stats_in, open(stats_out, "w") as stats_out:

    linecount = 0
    stats_dict = {}
    for line in stats_in:  
        if linecount == 0:
            linecount += 1
            stats_out.write(line.strip() + "\tperc_variant\tperc_perfect\tread_quality\ttelo_quality\n")
        else:
            stats_dict[line.strip().split()[0]] = line.strip()
    
    linecount = 0
    read = []
    for line in input_fh:
        linecount += 1
        read.append(line.strip())
        if linecount % 4 == 0:
            telo_seq = read[1][int(stats_dict[read[0].split()[0]].split()[3]):]
            telo_quality = get_mean_qual(read[3][int(stats_dict[read[0].split()[0]].split()[3]):])
            perfect_perc = telo_seq.count(repeat) * len(repeat) / len(telo_seq) * 100
            one_subs = len(list(regex.finditer(r'(%s){s<=1}' % repeat, telo_seq))) * len(repeat) / len(telo_seq) * 100
            stats_out.write('{}\t{}\t{}\t{}\t{}\n'.format(stats_dict[read[0].split()[0]], one_subs, perfect_perc, get_mean_qual(read[3]), telo_quality))
            read = []
        



    
