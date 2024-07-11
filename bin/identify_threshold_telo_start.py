import sys
import regex

input_file = sys.argv[1]
repeat = sys.argv[2]
perfect_repeats = sys.argv[3] == 'true'
sliding_window_size = int(sys.argv[4])
sliding_window_interval = int(sys.argv[5])
start_repeat_threshold = float(sys.argv[6])
end_repeat_threshold = float(sys.argv[7])
telomeric_repeat_percentage = float(sys.argv[8])
telo_reads_file = sys.argv[9]
non_telo_reads_file = sys.argv[10]

with open(telo_reads_file, 'w') as telo_fh, open(non_telo_reads_file, 'w') as non_telo_fh, open(input_file, 'r') as reads_fh:
    read = []
    linecount = 0 
    for line in reads_fh:
        linecount += 1
        if linecount % 4 == 0:
            read.append(line.strip())
            # determine telo start
            telo_found = False
            for i in range(0, len(read[1])-sliding_window_size, sliding_window_interval):
                subseq = read[1][i:i+sliding_window_size]
                if perfect_repeats:
                    telo_perc = (subseq.count(repeat) * len(repeat)) / sliding_window_size
                else:
                    telo_perc = len(list(regex.finditer(r"(%s){s<=1}" % repeat, subseq))) * len(repeat) / sliding_window_size
                
                if telo_perc >= start_repeat_threshold and not telo_found:
                    telo_found = True
                    if perfect_repeats:
                        telo_start = regex.search(repeat, subseq).span()[0] + i
                    else:
                        telo_start = regex.search(r'(%s){s<=1}' % repeat, subseq).span()[0] + i

                elif telo_perc >= end_repeat_threshold and telo_found:
                    pass
                else:
                    telo_found = False
            if telo_found and len(list(regex.finditer(r'(%s){s<=1}' % repeat, read[1][telo_start:]))) * len(repeat) / (len(read[1]) - telo_start) >= telomeric_repeat_percentage:
                    # write out to telo files
                telo_fh.write("{}\ttelo_start:{}\ttelo_length:{}\n{}\n{}\n{}\n".format(read[0], telo_start, len(read[1]) - telo_start, read[1], read[2], read[3]))
            else:
                    # write out to non telo file
                non_telo_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            read = []
        else:
            read.append(line.strip())
