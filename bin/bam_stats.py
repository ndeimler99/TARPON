#!/usr/bin/env python3

import argparse
import pysam
import numpy as np

def main(args):

    with open(args.out_file, "w") as out_fh:
        out_fh.write("file\tnum_seqs\tmean_length\tmin_length\tmax_length\tQ1\tQ2\tQ3\tmean_quality\tmin_quality\tmax_quality\tQ1_qual\tQ2_qual\tQ3_qual\n")

        for file in args.bam_files:
            print(file)
            aln_file = pysam.AlignmentFile(file, "rb", check_sq=False)
            
            num_seqs = 0
            seq_length = []
            seq_quality = []

            for aln in aln_file:
                num_seqs += 1
                seq_length.append(len(aln.query_sequence))
                seq_quality.append(np.mean(aln.query_qualities))
            
            if num_seqs > 0:

                out_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file, num_seqs, np.mean(seq_length), min(seq_length), max(seq_length), np.quantile(seq_length, 0.25),
                                                    np.quantile(seq_length, 0.5), np.quantile(seq_length, 0.75), np.mean(seq_quality), min(seq_quality), max(seq_quality),
                                                    np.quantile(seq_quality, 0.25), np.quantile(seq_quality, 0.5), np.quantile(seq_quality, 0.75)))
            else:
                out_fh.write("{}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n".format(file))
            aln_file.close()



def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_files", required=True, nargs="+")
    parser.add_argument("--out_file", required=True)
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
