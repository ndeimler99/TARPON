#!/usr/bin/env python3

import regex
import argparse
import pysam

def get_read_qual(seq):
    """Returns the mean phred converted quality score for a given seq"""
    return sum([i for i in seq]) / len(seq)

def get_telo_start(read, repeat, sliding_window, sliding_window_interval, upper_threshold, lower_threshold, consecutive_threshold, mutant=None):
    # identifies start of telomeric read (read).  See manuscript for more details on how this is done and why it is done each way
    telo_found = False
    telo_start = None
    below_threshold = 0
    for i in range(0, len(read)-sliding_window, sliding_window_interval):
        telo_matches = list(regex.finditer(r"(%s){s<=1}" % repeat, read[i:i+sliding_window]))
        if mutant is None:
            telo_perc = len(telo_matches) * len(repeat) / sliding_window
        else:
            telo_perc = (len(telo_matches) * len(repeat) + read[i:i+sliding_window].count(mutant) * len(mutant) ) / sliding_window

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

def check_valid(read, repeat, telomeric_rep_perc, mutant=None):
    # check if telomeric sequence from telo start to telo end is above telomeric_rep_perc
    
    if len(read) != 0 and read is not None:
        pass
    else:
        return False
    if mutant is None:
        telo_perc = len(list(regex.finditer(r"(%s){s<=1}" % repeat, read))) * len(repeat) / (len(read))
    else:
        telo_perc = (len(list(regex.finditer(r"(%s){s<=1}" % repeat, read))) * len(repeat) + read.count(mutant) * len(mutant)) / (len(read))

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
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--pre_telomeric_repeat_percentage", required=True)
    parser.add_argument("--pre_telo_distance", required=True)
    parser.add_argument("--minimum_telomere_length", required=True)

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
    args.pre_telomeric_repeat_percentage = float(args.pre_telomeric_repeat_percentage)
    args.pre_telo_distance = int(args.pre_telo_distance)
    args.minimum_telomere_length = int(args.minimum_telomere_length)

    input_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    telo_out = pysam.AlignmentFile(args.telomeric_fastq_out, "wb", template=input_fh)
    no_telo_out = pysam.AlignmentFile(args.no_telomere_out, "wb", template=input_fh)
    filtered_fh = pysam.AlignmentFile(args.filtered_out, "wb", template=input_fh)

    if args.mutant != "false":
        hamming = 0
        for i in range(0, len(args.mutant)):
            if args.mutant[i] != args.repeat[i]:
                hamming += 1
        if hamming <= 1:
            args.mutant = "false"
    
    # isolate through fastq file to perform telo start analysis on each individual read
    with open(args.stats_fh, "w") as stats_fh:
        stats_fh.write("read_id\tstrand\tread_len\tvrr_start_pos\tvrr_telo_length\ttelo_length\tread_qual\ttelo_qual\n")
        
        for aln in input_fh:
            # linecount += 1
            # if linecount % 4 == 0:
            #     read.append(line.strip())
            #     #analysis
            if args.mutant == "false":
                telo_found, telo_start = get_telo_start(aln.query_sequence, args.repeat, args.sliding_window, args.sliding_window_interval, args.upper_threshold, args.lower_threshold, args.consecutive_threshold)
            else:
                telo_found, telo_start = get_telo_start(aln.query_sequence, args.repeat, args.sliding_window, args.sliding_window_interval, args.upper_threshold, args.lower_threshold, args.consecutive_threshold, args.mutant)

            if not telo_found:
                    #write to file
                no_telo_out.write(aln)
                continue
            
            if args.mutant == "false":
                if not check_valid(aln.query_sequence[telo_start:], args.repeat, args.telomeric_rep_perc):
                    #write to file
                    filtered_fh.write(aln)
                    continue
                if telo_start - 2000 >= 0:
                    start_val = telo_start-2000
                else:
                    start_val = 0
                if check_valid(aln.query_sequence[start_val:telo_start], args.repeat, args.pre_telomeric_repeat_percentage):
                    filtered_fh.write(aln)
                    continue
            else:
                if not check_valid(aln.query_sequence[telo_start:], args.repeat, args.telomeric_rep_perc, args.mutant):
                    #write to file
                    filtered_fh.write(aln)
                    continue
                if telo_start - args.pre_telo_distance >= 0:
                    start_val = telo_start-args.pre_telo_distance
                else:
                    start_val = 0
                if check_valid(aln.query_sequence[start_val:telo_start], args.repeat, args.pre_telomeric_repeat_percentage, args.mutant):
                        filtered_fh.write(aln)
                        continue

            telo_length = get_telo_length(aln.query_sequence[telo_start:], args.repeat * args.consecutive_repeats)
            if len(aln.query_sequence) - telo_start < args.minimum_telomere_length:
                filtered_fh.write(aln)
                continue

            read_qual = get_read_qual(aln.query_qualities)
            telo_qual = get_read_qual(aln.query_qualities[telo_start:])
            #write to telo out fastq file
            telo_out.write(aln)
                #write to stats file
            stats_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(aln.query_name, aln.get_tag("XS"), len(aln.query_sequence), telo_start, len(aln.query_sequence)-telo_start, telo_length, read_qual, telo_qual))
            

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)