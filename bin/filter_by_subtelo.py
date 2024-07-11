#!/usr/bin/env python3

import sys
import regex

input_file = sys.argv[1]
min_subtelo_length = int(sys.argv[2])
min_subtelo_threshold = float(sys.argv[3])
repeat = sys.argv[4]
passes_subtelo = sys.argv[5]
fails_subtelo = sys.argv[6]


with open(input_file, 'r') as fh, \
 open(passes_subtelo, 'w') as subtelo_pass, open(fails_subtelo, 'w') as subtelo_fail:
    read = []
    linecount = 0
    for line in fh:
        linecount += 1
        if linecount % 4 == 0:
            read.append(line.strip())

            if len(read[1]) < min_subtelo_threshold:
                #write to subtelo fail
                subtelo_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            else:
                matches = list(regex.finditer(r'(%s){s<=1}' % repeat, read[1][0:min_subtelo_length]))
            #write to subtelo pass
                if len(matches) * len(repeat)/ min_subtelo_length <= min_subtelo_threshold:
                    subtelo_pass.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                else:
                    subtelo_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            read = []
        else:
            read.append(line.strip())
