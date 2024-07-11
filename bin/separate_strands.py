#!/usr/bin/env python3

import sys

input_file = sys.argv[1]
g_file = sys.argv[2]
c_file = sys.argv[3]

with open(input_file, 'r') as input_file_fh, open(g_file, 'w') as g_fh, open(c_file, 'w') as c_fh:
    linecount = 0
    read = []
    for line in input_file_fh:
        linecount += 1
        read.append(line.strip())
        if linecount % 4 == 0:
            if read[0].split()[1] == "G":
                g_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            else:
                c_fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            read = []
