import sys
import regex

input_file = sys.argv[1]
adaptor_sequence = sys.argv[2][0:12]
adaptor_errors = int(sys.argv[3])
out_file = sys.argv[4]
min_subtelo_length = int(sys.argv[5])
min_subtelo_threshold = float(sys.argv[6])
filtered_reads = sys.argv[7]
removed_reads = sys.argv[8]
repeat = sys.argv[9]

with open(input_file, 'r') as fh, open(out_file, 'w') as out_fh, open(filtered_reads, 'w') as filtered_fh, open(removed_reads, 'w') as removed_fh:
    read = []
    linecount = 0 
    for line in fh:
        linecount += 1
        if linecount % 4 == 0:
            read.append(line.strip())
            matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (adaptor_sequence, adaptor_errors), read[1]))
            if len(matches) > 0:
                for match in matches:
                    if read[1][0:match.span()[0]].count("GGTTAG") >= 10:
                        read[0] = read[0] + "\t" + read[1][match.span()[0]:match.span()[0] + 100]
                        read[1] = read[1][0:match.span()[0]]
                        read[3] = read[3][0:match.span()[0]]
                        out_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                        if len(read[1]) < min_subtelo_threshold:
                            removed_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                        else:
                            matches = list(regex.finditer(r'(%s){s<=1}' % repeat, read[1][0:min_subtelo_length]))
                            if len(matches) * len(repeat) / min_subtelo_length >= min_subtelo_threshold:
                                removed_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                            else:
                                filtered_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                        break  
            read = []
        else:
            read.append(line.strip())
