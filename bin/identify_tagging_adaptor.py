#!/usr/bin/env python3

import argparse
import regex
import pysam
import multiprocessing


def identify_adaptor(seq_dict, repeat, mutant, adaptor_sequence, adaptor_errors):
    
    adaptor_fail = []
    adaptor_out = {}
    for read in seq_dict:
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (adaptor_sequence, adaptor_errors), seq_dict[read]["seq"])) 
        if len(matches) >= 0:
            for match in matches:
                if mutant == "false":
                    repeat_count = seq_dict[read]["seq"][0:match.span()[0]].count(repeat)
                else:
                    repeat_count = seq_dict[read]["seq"][0:match.span()[0]].count(repeat) + seq_dict[read]["seq"][0:match.span()[0]].count(mutant)
                
                if repeat_count >= 20:
                    adaptor_out[read] = match.span()[0]
                    break
            else:
                adaptor_fail.append(read)
        else:
            adaptor_fail.append(read)

    return adaptor_fail, adaptor_out

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

def main(args):

    # convert parameters from strings to usable data types
    args.adaptor_errors = int(args.adaptor_errors)
    args.overhang_length = int(args.overhang_length)
    args.threads = int(args.threads)
    
    fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_out = pysam.AlignmentFile(args.adaptor_found, "wb", template=fh)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=fh)

    aln_dict = {}
    seq_dict = {}

    for aln in fh:
        aln_dict[aln.query_name] = aln
        seq_dict[aln.query_name] = {"seq":aln.query_sequence, "strand":aln.get_tag("XS")}
    
    # divide seq_dict into chunks
    bam_chunks = split_dict(seq_dict, args.threads)

    mp_args = [(chunk, args.repeat, args.mutant, args.adaptor_sequence, args.adaptor_errors) for chunk in bam_chunks]

    with multiprocessing.Pool(args.threads) as pool:
        results = pool.starmap(identify_adaptor, mp_args)

    for result in results:
        for read in result[0]:
            adaptor_fail.write(aln_dict[read])
        for read in result[1]:
            q = aln_dict[read].query_qualities
            if aln_dict[read].get_tag("XS") == "C":
                aln_dict[read].set_tag("XB", aln_dict[read].query_sequence[result[1][read]-args.overhang_length:result[1][read] + 100])
                aln_dict[read].query_sequence = aln_dict[read].query_sequence[0:result[1][read]-args.overhang_length]
                aln_dict[read].query_qualities = q[0:result[1][read]-args.overhang_length]
                adaptor_out.write(aln_dict[read])
            else:
                aln_dict[read].set_tag("XB", aln_dict[read].query_sequence[result[1][read]:result[1][read] + 100])
                aln_dict[read].query_sequence = aln_dict[read].query_sequence[0:result[1][read]]
                aln_dict[read].query_qualities = q[0:result[1][read]]
    # for aln in fh:
    #     matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.adaptor_sequence, args.adaptor_errors), aln.query_sequence)) 
    #     if len(matches) > 0:
    #         for match in matches:
    #                     # identify the first match that has at least ten telomeric repeats prior to the start of the match
    #                         # once this has been found chop the read at the start of the adaptor sequence and add the next 100bp to the header line and stop looping (discard portion of read further downstream)
                
    #             if args.mutant == "false":
    #                 repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) #/ aln.query_sequence.count(args.repeat) * 100
    #             else:
    #                 repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) + aln.query_sequence[0:match.span()[0]].count(args.mutant)

    #             if repeat_count >= 20:
    #                 if aln.get_tag("XS") == "C":
    #                     aln.set_tag("XB", aln.query_sequence[match.span()[0]-args.overhang_length:match.span()[0] + 100])
    #                 else:
    #                     aln.set_tag("XB", aln.query_sequence[match.span()[0]:match.span()[0] + 100])
    #                 q = aln.query_qualities
    #                 if aln.get_tag("XS") == "C":
    #                     #print(aln.query_sequence[match.span()[0]-50:match.span()[0]+50])
    #                     aln.query_sequence = aln.query_sequence[0:match.span()[0]-args.overhang_length]
    #                     #print(aln.query_sequence[match.span()[0]-200:match.span()[0]+50])
    #                     aln.query_qualities = q[0:match.span()[0]-args.overhang_length]
    #                     adaptor_out.write(aln)
    #                 else:
    #                     aln.query_sequence = aln.query_sequence[0:match.span()[0]]
    #                     aln.query_qualities = q[0:match.span()[0]]
    #                     adaptor_out.write(aln)
    #                 break
    #         else:
    #                     # this means no adaptor sequence was found after ten telomeric repeats
    #             print(aln.query_name)
    #             print(aln.query_sequence)
    #             adaptor_fail.write(aln)
    #     else:
    #                 # no adaptor found at all
    #         adaptor_fail.write(aln)
        
                
    fh.close()
    adaptor_out.close()
    adaptor_fail.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--adaptor_sequence", required=True)
    parser.add_argument("--adaptor_errors", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--adaptor_found", required=True)
    parser.add_argument("--no_adaptor", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--overhang_length", required=True)
    parser.add_argument("--threads", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
