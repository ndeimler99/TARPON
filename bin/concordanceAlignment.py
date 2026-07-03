#!/usr/bin/env python3

import argparse
from Bio import SeqIO, pairwise2

def main(args):

    sampleA = SeqIO.to_dict(SeqIO.parse(args.sampleA, "fasta"))
    sampleB = SeqIO.to_dict(SeqIO.parse(args.sampleB, "fasta"))

    sampleA_name = args.sampleA.split(".")[0]
    sampleB_name = args.sampleB.split(".")[0]
    
    with open("{}-{}.alignments.txt".format(sampleA_name, sampleB_name), "w") as out_fh:
        out_fh.write("{}\t{}\tp_ident\n".format(sampleA_name, sampleB_name))
        for readA in sampleA:
            for readB in sampleB:
                aln = pairwise2.align.globalxx(sampleA[readA].seq, sampleB[readB].seq, one_alignment_only=True)[0]
                matches = aln[2]
                length = max(len(sampleA[readA].seq), len(sampleB[readB].seq))
                percent_id = matches / length * 100
                out_fh.write("{}.{}\t{}.{}\t{}\n".format(sampleA_name, readA, sampleB_name, readB, percent_id))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sampleA", required=True)
    parser.add_argument("--sampleB", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
