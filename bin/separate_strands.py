#!/usr/bin/env python3

import argparse

def main(args):
    with open(args.input_file, 'r') as input_file_fh, open(args.g_file, 'w') as g_fh, open(args.c_file, 'w') as c_fh:
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

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--g_file", required=True)
    parser.add_argument("--c_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)