#!/usr/bin/env python3

import sys

sample_file = sys.argv[1]
barcode_errors = int(sys.argv[2])

barcodes = []
with open(sample_file, 'r') as sample_fh:
    linecount = 0
    for line in sample_fh:
        if linecount == 0:
            linecount += 1
            continue
        line = line.strip().split(",")
        barcodes.append(line[1])

for i in range(0, len(barcodes)-1):
    for j in range(i+1, len(barcodes)):
        if barcodes[i] == barcodes[j]:
            sys.exit(1)
        else:
            distance = 0
            for k in range(0, len(barcodes[i])):
                if barcodes[i][k] != barcodes[j][k]:
                    distance += 1
            if distance <= barcode_errors:
                sys.exit(1)

# f = open("passed.txt", "w")
# f.close()
            
        
        
