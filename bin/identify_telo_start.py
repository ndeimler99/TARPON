#!/usr/bin/env python3

import regex
import argparse

def get_read_qual(seq):
    """Returns the mean phred converted quality score for a given seq"""
    return sum([ord(i)-33 for i in seq]) / len(seq)

def get_telo_start(read, repeat, sliding_window, sliding_window_interval, upper_threshold, lower_threshold, consecutive_threshold):
    # identifies start of telomeric read (read).  See manuscript for more details on how this is done and why it is done each way
    telo_found = False
    telo_start = None
    below_threshold = 0
    for i in range(0, len(read)-sliding_window, sliding_window_interval):
        telo_matches = list(regex.finditer(r"(%s){s<=1}" % repeat, read[i:i+sliding_window]))
        telo_perc = len(telo_matches) * len(repeat) / sliding_window
        if telo_perc >= upper_threshold:
            below_threshold = 0
            if not telo_found:
                telo_found = True
                telo_start = i + telo_matches[0].span()[0]
        elif telo_perc < lower_threshold:
            if telo_found:
                below_threshold += 1
                if below_threshold >= consecutive_threshold:
                    telo_found = False
                    telo_start = None
    return telo_found, telo_start

def check_valid(read, repeat, telomeric_rep_perc):
    # check if telomeric sequence from telo start to telo end is above telomeric_rep_perc
    
    telo_perc = len(list(regex.finditer(r"(%s){s<=1}" % repeat, read))) * len(repeat) / (len(read))
    if telo_perc >= telomeric_rep_perc:
        return True
    else:
        return False

def get_telo_length(telomere, sequence):
    # return the telomere length which is currently is the number of nucleotides from telo start to telo end that are within 3 telomeric repeats in a row
    telo_matches = list(regex.finditer(sequence, telomere, overlapped=True))
    telomere = [*telomere]
    for match in telo_matches:
        for nucl in range(match.span()[0], match.span()[1]):
            telomere[nucl] = "N"
    return ''.join(telomere).count("N")

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--sliding_window", required=True, type=int)
    parser.add_argument("--sliding_window_interval", required=True, type=int)
    parser.add_argument("--upper_threshold", required=True, type=float)   
    parser.add_argument("--lower_threshold", required=True, type=float) 
    parser.add_argument("--telomeric_rep_perc", required=True, type=float)
    parser.add_argument("--consecutive_repeats", required=True, type=int)
    parser.add_argument("--consecutive_threshold", required=True, type=int)
    parser.add_argument("--telomeric_fastq_out", required=True)
    parser.add_argument("--no_telomere_out", required=True)
    parser.add_argument("--filtered_out", required=True)
    parser.add_argument("--stats_fh", required=True)

    
    return parser

def main(args):

    # convert parameters from strings to usable data types
    args.sliding_window = int(args.sliding_window)
    args.sliding_window_interval = int(args.sliding_window_interval)
    args.upper_threshold = float(args.upper_threshold)
    args.lower_threshold = float(args.lower_threshold)
    args.consecutive_repeats = int(args.consecutive_threshold)
    args.consecutive_repeats = int(args.consecutive_repeats)
    args.telomeric_rep_perc = float(args.telomeric_rep_perc)

    # isolate through fastq file to perform telo start analysis on each individual read
    with open(args.input_file, "r") as input_fh, open(args.telomeric_fastq_out, "w") as telo_out, \
    open(args.no_telomere_out, "w") as no_telo_out, open(args.filtered_out, 'w') as filtered_fh, \
    open(args.stats_fh, "w") as stats_fh:
        stats_fh.write("read_id\tstrand\tread_len\tvrr_start_pos\tvrr_telo_length\ttelo_length\tread_qual\ttelo_qual\n")
        read = []
        linecount = 0
        for line in input_fh:
            linecount += 1
            if linecount % 4 == 0:
                read.append(line.strip())
                #analysis
                telo_found, telo_start = get_telo_start(read[1], args.repeat, args.sliding_window, args.sliding_window_interval, args.upper_threshold, args.lower_threshold, args.consecutive_threshold)
                if not telo_found:
                    #write to file
                    no_telo_out.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                    read = []
                    continue
                
                if not check_valid(read[1][telo_start:], args.repeat, args.telomeric_rep_perc):
                    #write to file
                    filtered_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                    read = []
                    continue
                
                telo_length = get_telo_length(read[1][telo_start:], args.repeat * args.consecutive_repeats)
                read_qual = get_read_qual(read[3])
                telo_qual = get_read_qual(read[3][telo_start:])
                #write to telo out fastq file
                telo_out.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                #write to stats file
                stats_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read[0].split()[0], read[0].split()[1], len(read[1]), telo_start, len(read[1])-telo_start, telo_length, read_qual, telo_qual))
                read = []
            else:
                read.append(line.strip())

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)