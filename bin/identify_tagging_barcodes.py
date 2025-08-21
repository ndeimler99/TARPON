#!/usr/bin/env python3

import argparse
import regex
import pysam 
import multiprocessing

def identify_first_barcode(read, barcode_dict, barcode_errors, repeat, mutant):
    """identifies location of the first barcode existing within a given read after ten telomeric repeats are identified:
    returns none if no barcodes are found. returns the first barcode and its location when found"""
    location_dict = {}
    for sample in barcode_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (barcode_dict[sample], barcode_errors), read))
        if len(matches) > 0:
            for match in matches:
                if mutant == "false":
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

def split_dict(d, n):
        """Split dictionary d into n smaller dictionaries."""
        items = list(d.items())
        k, m = divmod(len(items), n)
        chunks = []
        for i in range(n):
            start = i*k + min(i, m)
            end = (i+1)*k + min(i+1, m)
            chunks.append(dict(items[start:end]))
        return chunks

def identify_barcodes(seq_dict, barcode_dict, barcode_errors, repeat, mutant):

    adaptor_fail = []
    adaptor_loc = {}
    sample_dict = {}

    for read in seq_dict:
        barcode, location = identify_first_barcode(seq_dict[read], barcode_dict, barcode_errors, repeat, mutant)
        if barcode is None:
            adaptor_fail.append(read)
            continue
        else:
            adaptor_loc[read] = location
            sample_dict[read] = barcode
    
    return adaptor_fail, adaptor_loc, sample_dict


def main(args):

    args.barcode_errors = int(args.barcode_errors)
    args.overhang_length = int(args.overhang_length)
    args.threads = int(args.threads)

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
    
    aln_dict = {}
    seq_dict = {}

    for aln in in_fh:
        aln_dict[aln.query_name] = aln
        seq_dict[aln.query_name] = aln.query_sequence
    
    # divide seq_dict into chunks
    bam_chunks = split_dict(seq_dict, args.threads)

    mp_args = [(chunk, barcode_dict, args.barcode_errors, args.repeat, args.mutant) for chunk in bam_chunks]

    with multiprocessing.Pool(args.threads) as pool:
        results = pool.starmap(identify_barcodes, mp_args)

    for result in results:
        for read in result[0]:
            adaptor_fail.write(aln_dict[read])
        for read in result[1]:
            # this is the adaptor location
            q = aln_dict[read].query_qualities
            if aln_dict[read].get_tag("XS") == "C":
                aln_dict[read].set_tag("XB", aln_dict[read].query_sequence[result[1][read]-args.overhang_length:result[1][read] + 100])
                aln_dict[read].query_sequence = aln_dict[read].query_sequence[0:result[1][read]-args.overhang_length]
                aln_dict[read].query_qualities = q[0:result[1][read]-args.overhang_length]
            else:
                aln_dict[read].set_tag("XB", aln_dict[read].query_sequence[result[1][read]:result[1][read] + 100])
                aln_dict[read].query_sequence = aln_dict[read].query_sequence[0:result[1][read]]
                aln_dict[read].query_qualities = q[0:result[1][read]]
            read_dict[result[2][read]].append(read)
    # #with open(args.input_file, 'r') as fh, open(args.no_adaptor, 'w') as adaptor_fail:
    #  #   read = []
    #  #   linecount = 0 
    # for aln in in_fh:
    #     #linecount += 1
    #     # iterate through reads in a fastq file. Identify barcode for every read
    #     # add barcode sequence to header row for documentation
    #     #if linecount % 4 == 0:
    #         #read.append(line.strip())
        
    #     if args.mutant == "false":
    #         barcode, location = identify_first_barcode(aln.query_sequence, barcode_dict, args.barcode_errors, args.repeat)
    #     else:
    #         print("within else statement, mutant was specified")
    #         barcode, location = identify_first_barcode(aln.query_sequence, barcode_dict, args.barcode_errors, args.repeat, args.mutant)

    #     if barcode is not None and location != len(aln.query_sequence):
    #         q = aln.query_qualities
    #         if aln.get_tag("XS") == "C":
    #             aln.set_tag("XB", aln.query_sequence[location-args.overhang_length:location+100])
    #             aln.query_sequence = aln.query_sequence[0:location-args.overhang_length]
    #             aln.query_qualities = q[0:location-args.overhang_length]
    #         else:
    #             aln.set_tag("XB", aln.query_sequence[location:location+100])
    #             aln.query_sequence = aln.query_sequence[0:location]
    #             aln.query_qualities = q[0:location]
    #         read_dict[barcode].append(aln)
    #     else:
    #         adaptor_fail.write(aln)

    # write out each individual fastq file from demultiplexing
    for sample in read_dict:
        fh = pysam.AlignmentFile("{}/{}.bam".format(args.out_fh, sample), "wb", template=in_fh)
        for read in read_dict[sample]:
            fh.write(aln_dict[read])
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
    parser.add_argument("--threads", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
