#!/usr/bin/env python3

import argparse
import regex
import pysam

def identify_first_barcode(header, barcode_dict, barcode_errors, repeat):
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
    args.adaptor_errors = int(args.adaptor_errors)
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

    # iterate through input_file by read checking header line for adaptor sequence
    # identify the first barcode in the header sequence
    # return its identity
    in_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=in_fh)

    for aln in in_fh:

        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.adaptor_sequence, args.adaptor_errors), aln.query_sequence)) 
        if len(matches) > 0:
            for match in matches:
                        # identify the first match that has at least ten telomeric repeats prior to the start of the match
                            # once this has been found chop the read at the start of the adaptor sequence and add the next 100bp to the header line and stop looping (discard portion of read further downstream)
                if args.mutant == "false":
                    repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat)
                else:
                    repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) + aln.query_sequence[0:match.span()[0]].count(args.mutant)
                if repeat_count >= 20:
                    q = aln.query_qualities

                    if aln.get_tag("XS") == "C":
                        barcode_seq = aln.query_sequence[match.span()[0]-args.overhang_length:match.span()[0] + 100]
                        aln.set_tag("XB", barcode_seq)
                        aln.query_sequence = aln.query_sequence[0:match.span()[0]-args.overhang_length]
                        aln.query_qualities = q[0:match.span()[0]-args.overhang_length]
                    else:
                        barcode_seq = aln.query_sequence[match.span()[0]:match.span()[0] + 100]
                        aln.set_tag("XB", barcode_seq)
                        aln.query_sequence = aln.query_sequence[0:match.span()[0]]
                        aln.query_qualities = q[0:match.span()[0]]

                    
                    barcode, location = identify_first_barcode(barcode_seq, barcode_dict, args.barcode_errors, args.repeat)
                    if barcode is not None and location[0] != len(barcode_seq):
                        barcode_seq = barcode_seq[location[0]:location[1]]
                        aln.set_tag("XB", barcode_seq)
                        read_dict[barcode].append(aln)
                        break
                    else:
                        # no barcode found
                        adaptor_fail.write(aln)
            else:
                            # this means no adaptor sequence was found after ten telomeric repeats
                adaptor_fail.write(aln)
        else:
                    # no adaptor found at all
            adaptor_fail.write(aln)

    # if it does not exist append to no adaptor file
    # else add to second dictionary

    # write out to sample files
    for sample in read_dict:
        fh = pysam.AlignmentFile("{}/{}.bam".format(args.out_prefix, sample), "wb", template=in_fh)
        for read in read_dict[sample]:
            fh.write(aln)
        fh.close()
    
    in_fh.close()
    adaptor_fail.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--sample_file", required=True)
    parser.add_argument("--barcode_errors", required=True)
    parser.add_argument("--no_adaptor", required=True)
    parser.add_argument("--out_prefix", required=True)
    parser.add_argument("--adaptor_sequence", required=True)
    parser.add_argument("--adaptor_errors", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--overhang_length"), required=True
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)