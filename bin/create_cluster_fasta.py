#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd

def main(args):

    # load cluster dict from out_file
    stats_df = pd.read_table(args.telo_stats, sep="\t")
    aln_file = pysam.AlignmentFile(args.telobam, "r", check_sq=False)

    aln_dict = {}
    for aln in aln_file:
        aln_dict[aln.query_name] = aln

    for cluster in stats_df["Cluster"].dropna().unique():
        with open("cluster_{}.fa".format(int(cluster)), "w") as fasta_fh:
            for read in stats_df.loc[stats_df["Cluster"]==cluster]["read_id"]:
                telo_start = stats_df.loc[stats_df["read_id"]==read, "vrr_start_pos"].iloc[0]
                if telo_start - 1000 < 0:
                    fasta_fh.write(">{}\n{}\n".format(read, aln_dict[read].query_sequence[0:telo_start + 2000]))
                else:
                    fasta_fh.write(">{}\n{}\n".format(read, aln_dict[read].query_sequence[telo_start-1000:telo_start+2000]))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telo_stats",required=True)
    parser.add_argument("--telobam", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)