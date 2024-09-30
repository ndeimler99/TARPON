#!/usr/bin/env python3

import argparse
import regex

def main(args):

    args.adaptor_errors = int(args.adaptor_errors)
    
    with open(args.input_file, 'r') as fh, \
    open(args.adaptor_found, 'w') as adaptor_out, open(args.no_adaptor, 'w') as adaptor_fail:
        read = []
        linecount = 0 
        for line in fh:
            linecount += 1
            if linecount % 4 == 0:
                read.append(line.strip())
                
                matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.adaptor_sequence, args.adaptor_errors), read[1])) 
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

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--adaptor_sequence", required=True)
    parser.add_argument("--adaptor_errors", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--adaptor_found", required=True)
    parser.add_argument("--no_adaptor", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
