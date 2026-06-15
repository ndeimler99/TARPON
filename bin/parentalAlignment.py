#!/usr/bin/env python3

import argparse
from Bio import SeqIO, pairwise2

def main(args):

    offspring_fa_dict = SeqIO.to_dict(SeqIO.parse(args.offspring, "fasta"))
    parental_fa_dict = SeqIO.to_dict(SeqIO.parse(args.parental, "fasta"))

    with open(args.output, "w") as out_fh:
        out_fh.write("offspring_cluster\tparental_cluster\tp_ident\n")
        for offspring_read in offspring_fa_dict:
            for parental_read in parental_fa_dict:
                aln = pairwise2.align.globalxx(offspring_fa_dict[offspring_read].seq, parental_fa_dict[parental_read].seq, one_alignment_only=True)[0]
                matches = aln[2]
                length = max(len(offspring_fa_dict[offspring_read].seq), len(parental_fa_dict[parental_read].seq))
                percent_id = matches / length * 100
                out_fh.write("{}\t{}\t{}\n".format(offspring_read, parental_read, percent_id))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--offspring", required=True)
    parser.add_argument("--parental", required=True)
    parser.add_argument("--output", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
