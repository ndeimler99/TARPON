#!/usr/bin/env python3

import argparse
import regex
import pysam

def main(args):

    # convert parameters from strings to usable data types
    args.adaptor_errors = int(args.adaptor_errors)
    args.overhang_length = int(args.overhang_length)
    
    fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_out = pysam.AlignmentFile(args.adaptor_found, "wb", template=fh)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=fh)
                
    for aln in fh:
    
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.adaptor_sequence, args.adaptor_errors), aln.query_sequence)) 
        if len(matches) > 0:
            for match in matches:
                        # identify the first match that has at least ten telomeric repeats prior to the start of the match
                            # once this has been found chop the read at the start of the adaptor sequence and add the next 100bp to the header line and stop looping (discard portion of read further downstream)
                if args.mutant == "false":
                    repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) / aln.query_sequence.count(args.repeat) * 100
                else:
                    repeat_count = aln.query_sequence[0:match.span()[0]].count(args.repeat) + aln.query_sequence[0:match.span()[0]].count(args.mutant)

                if repeat_count >= 20:
                    if aln.get_tag("XS") == "C":
                        aln.set_tag("XB", aln.query_sequence[match.span()[0]-args.overhang_length:match.span()[0] + 100])
                    else:
                        aln.set_tag("XB", aln.query_sequence[match.span()[0]:match.span()[0] + 100])
                    q = aln.query_qualities
                    if aln.get_tag("XS") == "C":
                        #print(aln.query_sequence[match.span()[0]-50:match.span()[0]+50])
                        aln.query_sequence = aln.query_sequence[0:match.span()[0]-args.overhang_length]
                        print(aln.query_sequence[match.span()[0]-200:match.span()[0]+50])
                        aln.query_qualities = q[0:match.span()[0]-args.overhang_length]
                        adaptor_out.write(aln)
                    else:
                        aln.query_sequence = aln.query_sequence[0:match.span()[0]]
                        aln.query_qualities = q[0:match.span()[0]]
                        adaptor_out.write(aln)
                    break
            else:
                        # this means no adaptor sequence was found after ten telomeric repeats
                adaptor_fail.write(aln)
        else:
                    # no adaptor found at all
            adaptor_fail.write(aln)
        
                
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
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
