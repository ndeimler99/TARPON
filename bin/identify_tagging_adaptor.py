#!/usr/bin/env python3

import sys
import regex

input_file = sys.argv[1]
adaptor_sequence = sys.argv[2]
adaptor_errors = int(sys.argv[3])
repeat = sys.argv[4]
adaptor_found = sys.argv[5]
no_adaptor = sys.argv[6]

with open(input_file, 'r') as fh, \
 open(adaptor_found, 'w') as adaptor_out, open(no_adaptor, 'w') as adaptor_fail:
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
                        adaptor_out.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                        break
                else:
                    adaptor_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
            else:
                # write to adaptor not found
                adaptor_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
    
            read = []
        else:
            read.append(line.strip())
