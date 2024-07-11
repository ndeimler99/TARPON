#!/usr/bin/env python3

import sys
import regex

input_file = sys.argv[1]
repeat = sys.argv[2]
sliding_window = int(sys.argv[3])
sliding_window_interval = int(sys.argv[4])
upper_threshold = float(sys.argv[5])
lower_threshold = float(sys.argv[6])
telomeric_rep_perc = float(sys.argv[7])
consecutive_repeats = int(sys.argv[8])
telomeric_fastq_out = sys.argv[9]
no_telomere_out = sys.argv[10]
filtered_out = sys.argv[11]
stats_fh = sys.argv[12]


def get_telo_start(read, repeat, sliding_window, sliding_window_interval, upper_threshold, lower_threshold):
    telo_found = False
    telo_start = None
    for i in range(0, len(read)-sliding_window, sliding_window_interval):
        telo_matches = list(regex.finditer(r"(%s){s<=1}" % repeat, read[i:i+sliding_window]))
        telo_perc = len(telo_matches) * len(repeat) / sliding_window
        if telo_perc >= upper_threshold:
            if not telo_found:
                telo_found = True
                telo_start = i + telo_matches[0].span()[0]
        elif telo_perc < lower_threshold:
            if telo_found:
                telo_found = False
                telo_start = None
    
    return telo_found, telo_start

def check_valid(read, repeat, telomeric_rep_perc):
    
    telo_perc = len(list(regex.finditer(r"(%s){s<=1}" % repeat, read))) * len(repeat) / (len(read))
    if telo_perc >= telomeric_rep_perc:
        return True
    else:
        return False

def get_telo_length(telomere, sequence):
    telo_matches = list(regex.finditer(sequence, telomere, overlapped=True))
    telomere = [*telomere]
    for match in telo_matches:
        for nucl in range(match.span()[0], match.span()[1]):
            telomere[nucl] = "N"
    return ''.join(telomere).count("N")


with open(input_file, "r") as input_fh, open(telomeric_fastq_out, "w") as telo_out, \
 open(no_telomere_out, "w") as no_telo_out, open(filtered_out, 'w') as filtered_fh, \
 open(stats_fh, "w") as stats_fh:
    stats_fh.write("read_id\tstrand\tread_len\tvrr_start_pos\tvrr_telo_length\ttelo_length\n")
    read = []
    linecount = 0
    for line in input_fh:
        linecount += 1
        if linecount % 4 == 0:
            read.append(line.strip())
            #analysis
            telo_found, telo_start = get_telo_start(read[1], repeat, sliding_window, sliding_window_interval, upper_threshold, lower_threshold)
            if not telo_found:
                #write to file
                no_telo_out.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                read = []
                continue
            
            if not check_valid(read[1][telo_start:], repeat, telomeric_rep_perc):
                #write to file
                filtered_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                read = []
                continue
            
            telo_length = get_telo_length(read[1][telo_start:], repeat * consecutive_repeats)
            #write to telo out fastq file
            telo_out.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            #write to stats file
            stats_fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read[0].split()[0], read[0].split()[1], len(read[1]), telo_start, len(read[1])-telo_start, telo_length))
            read = []
        else:
            read.append(line.strip())
