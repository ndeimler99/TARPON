#!/usr/bin/env python3

import argparse
import regex
import pysam

def main(args):

    # convert parameters from strings to usable data types    
    fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    adaptor_out = pysam.AlignmentFile(args.adaptor_found, "wb", template=fh)
    adaptor_fail = pysam.AlignmentFile(args.no_adaptor, "wb", template=fh)
                
    for aln in fh:
    
        matches = list(regex.finditer(r'(?e)(%s){e<=%s}' % (args.repeat*2, 2), aln.query_sequence)) 
        if len(matches) > 0:
            last_match = matches[-1]
            aln.set_tag("XB", aln.query_sequence[last_match.span()[0]:])
            q = aln.query_qualities
            aln.query_sequence = aln.query_sequence[0:last_match.span()[0]]
            aln.query_qualities = q[0:last_match.span()[0]]
            adaptor_out.write(aln)
        else:
                    # no adaptor found at all
            adaptor_fail.write(aln)
        
    fh.close()
    adaptor_out.close()
    adaptor_fail.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--adaptor_found", required=True)
    parser.add_argument("--no_adaptor", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
