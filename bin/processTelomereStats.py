#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat_files", nargs="+", required=True)

    #    parser.add_argument("--restriction_digest_analysis", required=False, default="test")
    return parser

def main(args):
    # convert parameters from strings to usable data types

    vrr_dict = {}
    out_df = []

    # open stats files
    vrr_fh = open("combined_stats.VRR.txt", "w")

    vrr_fh.write("Sample_ID\tNumber_of_Reads\tMean_VRR_Telomere_Length\tStandard_Deviation_VRR_Length\tQ1\tQ2\tQ3\tMin_VRR_Telo_Length\tMax_VRR_Telo_Length\tmin_qual\tmax_qual\tqual_q1\tqual_q2\tqual_q3\tmean_qual\n")

    # for every sample in demultiplexed data write out telomere statistics
    for file in args.stat_files:
        df = pd.read_table(file, sep="\t")
        # stats
        df["sample"] = file.strip().split("/")[-1].split(".")[0]

        vrr_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['vrr_telo_length']), \
                                                                np.mean(df['vrr_telo_length']), np.std(df['vrr_telo_length']), np.quantile(df['vrr_telo_length'], 0.25), np.quantile(df['vrr_telo_length'], 0.5), \
                                                                    np.quantile(df['vrr_telo_length'], 0.75), min(df['vrr_telo_length']), max(df['vrr_telo_length']), min(df['read_qual']), max(df['read_qual']), \
                                                                    np.quantile(df['read_qual'], 0.25), np.quantile(df["read_qual"], 0.5), np.quantile(df["read_qual"], 0.75), np.mean(df["read_qual"])))

    vrr_fh.close()

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
