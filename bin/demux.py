#!/usr/bin/env python3

import argparse
import regex

def identify_first_barcode(header, barcode_dict, barcode_errors):
    """identifies location of the first barcode existing within a given read after ten telomeric repeats are identified:
    returns none if no barcodes are found. returns the first barcode and its location when found"""
    location_dict = {}
    for sample in barcode_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (barcode_dict[sample], barcode_errors), header))
        if len(matches) > 0:
            location_dict[sample] = [matches[0].span()[0], matches[0].span()[1]]
    
    if len(location_dict) == 0:
        return None, None
    
    min_loc = [len(header), len(header)]
    barcode = None
    for sample in location_dict:
        if location_dict[sample][0] < min_loc[0]:
            min_loc = location_dict[sample]
            barcode = sample
    return barcode, min_loc


def main(args):

    args.barcode_errors = int(args.barcode_errors)

    barcode_dict = {}
    read_dict = {}
    first_line = True
    # create a barcode dictionary
    with open(args.sample_file, 'r') as sample_file:
        for line in sample_file:
            if first_line:
                first_line = False
                continue
            line = line.strip().split(",")
            barcode_dict[line[0]] = line[1]
            read_dict[line[0]] = []

    # iterate through input_file by read checking header line for adaptor sequence
    # identify the first barcode in the header sequence
    # return its identity
    with open(args.input_file, 'r') as fh, open(args.no_adaptor, 'a') as adaptor_fail:
        read = []
        linecount = 0 
        for line in fh:
            linecount += 1
            # iterate through reads in a fastq file. Identify barcode for every read
            # add barcode sequence to header row for documentation
            if linecount % 4 == 0:
                read.append(line.strip())
                header = read[0].split("\t")[-1]
                print(read[0])
                print(header)
                barcode, location = identify_first_barcode(header, barcode_dict, args.barcode_errors)
                if barcode is not None and location[0] != len(header):
                    read_dict[barcode].append(read)
                    read[0] = read[0].split("\t")[:-1]
                    read[0].append(header[location[0]:location[1]])
                    read[0] = read[0].join("\t")
                    print(read[0])
                else:
                    adaptor_fail.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))
                read = []
            else:
                read.append(line.strip())

    # if it does not exist append to no adaptor file
    # else add to second dictionary

    # write out to sample files
    for sample in read_dict:
        with open("{}/{}.fastq".format(args.out_prefix, sample), 'w') as fh:
            for read in read_dict[sample]:
                fh.write("{}\n{}\n{}\n{}\n".format(read[0], read[1], read[2], read[3]))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--barcode_errors", required=True)
    parser.add_argument("--no_adaptor", required=True)
    parser.add_argument("--out_prefix", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)