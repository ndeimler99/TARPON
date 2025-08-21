#!/usr/bin/env python3

import argparse
import regex
import pysam
import multiprocessing

def identify_first_barcode(read, adaptor_start, barcode_dict, barcode_errors):

    """identifies location of the first barcode existing within a given read after ten telomeric repeats are identified:
    returns none if no barcodes are found. returns the first barcode and its location when found"""
    location_dict = {}
    for sample in barcode_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (barcode_dict[sample], barcode_errors), read))
        if len(matches) > 0:
            location_dict[sample] = [matches[0].span()[0], matches[0].span()[1]]
    
    if len(location_dict) == 0:
        return None, None
    
    min_loc = [len(read), len(read)]
    barcode = None
    for sample in location_dict:
        if location_dict[sample][0] < min_loc[0]:
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


def identify_adaptor(seq, repeat, mutant, adaptor_sequence, adaptor_errors):
    
    adaptor_start = None

    matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (adaptor_sequence, adaptor_errors), seq)) 
    if len(matches) >= 0:
        for match in matches:
            if mutant == "false":
                repeat_count = seq[0:match.span()[0]].count(repeat)
            else:
                repeat_count = seq[0:match.span()[0]].count(repeat) + seq[0:match.span()[0]].count(mutant)
                
            if repeat_count >= 20:
                adaptor_start = match.span()[0]
                break

    return adaptor_start

def identify_adaptor_and_demux(seq_dict, repeat, mutant, adaptor_sequence, adaptor_errors, barcode_dict, barcode_errors):
    
    adaptor_fail = []
    adaptor_loc = {}
    sample_dict = {}

    for read in seq_dict:
        adaptor_start = identify_adaptor(seq_dict[read], repeat, mutant, adaptor_sequence, adaptor_errors)

        if adaptor_start is None:
            adaptor_fail.append(read)
            continue
        else:
            barcode, barcode_pos = identify_first_barcode(seq_dict[read][adaptor_start:adaptor_start+200], adaptor_start, barcode_dict, barcode_errors)
            if barcode is None:
                adaptor_fail.append(read)
                continue
            else:
                sample_dict[read] = barcode
                adaptor_loc[read] = adaptor_start

    return adaptor_fail, adaptor_loc, sample_dict

def main(args):

    args.barcode_errors = int(args.barcode_errors)
    args.adaptor_errors = int(args.adaptor_errors)
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

    # iterate through input_file by read checking header line for adaptor sequence
    # identify the first barcode in the header sequence
    # return its identity
    in_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=in_fh)

    aln_dict = {}
    seq_dict = {}

    for aln in in_fh:
        aln_dict[aln.query_name] = aln
        seq_dict[aln.query_name] = aln.query_sequence
    
    # divide seq_dict into chunks
    bam_chunks = split_dict(seq_dict, args.threads)

    mp_args = [(chunk, args.repeat, args.mutant, args.adaptor_sequence, args.adaptor_errors, barcode_dict, args.barcode_errors) for chunk in bam_chunks]

    with multiprocessing.Pool(args.threads) as pool:
        results = pool.starmap(identify_adaptor_and_demux, mp_args)

    for result in results:
        for read in result[0]:
            adaptor_fail.write(aln_dict[read])
        for read in result[1]:
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

    # for aln in in_fh:
    #     print(aln.query_name)
    #     matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.adaptor_sequence, args.adaptor_errors), aln.query_sequence)) 
    #     if len(matches) > 0:
    #         for match in matches:
    #                     # identify the first match that has at least ten telomeric repeats prior to the start of the match
    #                         # once this has been found chop the read at the start of the adaptor sequence and add the next 100bp to the header line and stop looping (discard portion of read further downstream)
    #             if args.mutant == "false":
    #                 repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat)
    #             else:
    #                 repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) + aln.query_sequence[0:match.span()[0]].count(args.mutant)
    #             if repeat_count >= 20:
    #                 q = aln.query_qualities

    #                 if aln.get_tag("XS") == "C":
    #                     barcode_seq = aln.query_sequence[match.span()[0]-args.overhang_length:match.span()[0] + 200]
    #                     aln.set_tag("XB", barcode_seq)
    #                     aln.query_sequence = aln.query_sequence[0:match.span()[0]-args.overhang_length]
    #                     aln.query_qualities = q[0:match.span()[0]-args.overhang_length]
    #                 else:
    #                     barcode_seq = aln.query_sequence[match.span()[0]:match.span()[0] + 200]
    #                     aln.set_tag("XB", barcode_seq)
    #                     aln.query_sequence = aln.query_sequence[0:match.span()[0]]
    #                     aln.query_qualities = q[0:match.span()[0]]

                    
    #                 barcode, location = identify_first_barcode(barcode_seq, barcode_dict, args.barcode_errors)
    #                 if barcode is not None and location[0] != len(barcode_seq):
    #                     barcode_seq = barcode_seq[location[0]:location[1]]
    #                     aln.set_tag("XB", barcode_seq)
    #                     read_dict[barcode].append(aln)
    #                     break
    #                 else:
    #                     # no barcode found
    #                     adaptor_fail.write(aln)
    #         else:
    #                         # this means no adaptor sequence was found after ten telomeric repeats
    #             adaptor_fail.write(aln)
    #     else:
    #                 # no adaptor found at all
    #         adaptor_fail.write(aln)

    # if it does not exist append to no adaptor file
    # else add to second dictionary

    # write out to sample files
    for sample in read_dict:
        fh = pysam.AlignmentFile("{}/{}.bam".format(args.out_prefix, sample), "wb", template=in_fh)
        for read in read_dict[sample]:
            fh.write(aln_dict[read])
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
    parser.add_argument("--overhang_length", required=True)
    parser.add_argument('--threads', required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
