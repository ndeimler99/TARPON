import sys
import regex

input_file = sys.argv[1]
repeat = sys.argv[2]
repeat_sequence = int(sys.argv[3]) * repeat
telomeric_repeat_percentage = float(sys.argv[4])
telo_reads_file = sys.argv[5]
non_telo_reads_file = sys.argv[6]

with open(telo_reads_file, 'w') as telo_fh, open(non_telo_reads_file, 'w') as non_telo_fh, open(input_file, 'r') as reads_fh:
    read = []
    linecount = 0 
    for line in reads_fh:
        linecount += 1
        if linecount % 4 == 0:
            read.append(line.strip())
            # determine telo start
            telo_found = False
            match = regex.search(read[1], repeat_sequence)
            if match is not None:
                telo_found = True
                telo_start = match.span()[0]
            if telo_found and len(list(regex.finditer(r'(%s){s<=1}' % repeat, read[1][telo_start:]))) * len(repeat) / (len(read[1]) - telo_start) >= telomeric_repeat_percentage:
                    # write out to telo files
                telo_fh.write("{}\ttelo_start:{}\ttelo_length:{}\n{}\n{}\n{}\n".format(read[0], telo_start, len(read[1]) - telo_start, read[1], read[2], read[3]))
            else:
                    # write out to non telo file
                non_telo_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            read = []
        else:
            read.append(line.strip())
