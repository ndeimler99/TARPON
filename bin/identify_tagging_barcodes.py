#!/usr/bin/env python3

import argparse
import regex
import pysam 

def identify_first_barcode(read, barcode_dict, barcode_errors, repeat, mutant=None):
    """identifies location of the first barcode existing within a given read after ten telomeric repeats are identified:
    returns none if no barcodes are found. returns the first barcode and its location when found"""
    location_dict = {}
    for sample in barcode_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (barcode_dict[sample], barcode_errors), read))
        if len(matches) > 0:
            for match in matches:
                if mutant is None:
                    repeat_count = read[0:match.span()[0]].count(repeat)
                else:
                    repeat_count = read[0:match.span()[0]].count(repeat) + read[0:match.span()[0]].count(mutant)

                if repeat_count >= 20:
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
    args.overhang_length = int(args.overhang_length)

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

    in_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=in_fh)

    #with open(args.input_file, 'r') as fh, open(args.no_adaptor, 'w') as adaptor_fail:
     #   read = []
     #   linecount = 0 
    for aln in in_fh:
        #linecount += 1
        # iterate through reads in a fastq file. Identify barcode for every read
        # add barcode sequence to header row for documentation
        #if linecount % 4 == 0:
            #read.append(line.strip())
        
        if args.mutant == "false":
            barcode, location = identify_first_barcode(aln.query_sequence, barcode_dict, args.barcode_errors, args.repeat)
        else:
            print("within else statement, mutant was specified")
            barcode, location = identify_first_barcode(aln.query_sequence, barcode_dict, args.barcode_errors, args.repeat, args.mutant)

        if barcode is not None and location != len(aln.query_sequence):
            q = aln.query_qualities
            if aln.get_tag("XS") == "C":
                aln.query_sequence = aln.query_sequence[0:location-args.overhang_length]
                aln.query_qualities = q[0:location-args.overhang_length]
                aln.set_tag("XB", aln.query_sequence[location-args.overhang_length:location+100])
            else:
                aln.query_sequence = aln.query_sequence[0:location]
                aln.query_qualities = q[0:location]
                aln.set_tag("XB", aln.query_sequence[location:location+100])
            read_dict[barcode].append(aln)
        else:
            adaptor_fail.write(aln)

    # write out each individual fastq file from demultiplexing
    for sample in read_dict:
        fh = pysam.AlignmentFile("{}/{}.bam".format(args.out_fh, sample), "wb", template=in_fh)
        for read in read_dict[sample]:
            fh.write(read)
        fh.close()

    adaptor_fail.close()
    in_fh.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--barcode_errors", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--out_fh", required=True)
    parser.add_argument("--no_adaptor", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--overhang_length", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
