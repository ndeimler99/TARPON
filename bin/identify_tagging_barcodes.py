#!/usr/bin/env python3

import argparse
import regex

def identify_first_barcode(read, barcode_dict, barcode_errors, repeat):
    location_dict = {}
    for sample in barcode_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (barcode_dict[sample], barcode_errors), read))
        if len(matches) > 0:
            for match in matches:
                if read[0:match.span()[0]].count(repeat) >= 10:
                    location_dict[sample] = match.span()[0]
    
    if len(location_dict) == 0:
        return None, None
    
    min_loc = len(read)
    barcode = None
    for sample in location_dict:
        if location_dict[sample] < min_loc:
            min_loc = location_dict[sample]
            barcode = sample
    return barcode, min_loc

def main(args):

    args.barcode_errors = int(args.barcode_errors)

    barcode_dict = {}
    read_dict = {}
    first_line = True
    with open(args.sample_file, 'r') as sample_file:
        for line in sample_file:
            if first_line:
                first_line = False
                continue
            line = line.strip().split(",")
            barcode_dict[line[0]] = line[1]
            read_dict[line[0]] = []

    with open(args.input_file, 'r') as fh, open(args.no_adaptor, 'w') as adaptor_fail:
        read = []
        linecount = 0 
        for line in fh:
            linecount += 1
            if linecount % 4 == 0:
                read.append(line.strip())
                
                barcode, location = identify_first_barcode(read[1], barcode_dict, args.barcode_errors, args.repeat)
                if barcode is not None and location != len(read[1]):
                    #shorten read and quality score and write out to appropriate barcode file
                    #include 100 bp downstream in header line
                    read[1] = read[1][0:location]
                    read[3] = read[3][0:location]
                    read[0] = read[0] + "\t" + read[1][location:location+100]
                    read_dict[barcode].append(read)
                else:
                    adaptor_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                read = []
            else:
                read.append(line.strip())

    for sample in read_dict:
        with open("{}/{}.fastq".format(args.out_fh, sample), 'w') as fh:
            for read in read_dict[sample]:
                fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--barcode_errors", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--out_fh", required=True)
    parser.add_argument("--no_adaptor", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
