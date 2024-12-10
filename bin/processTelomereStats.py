#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat_files", nargs="+", required=True)
    parser.add_argument("--vrr_length", required=True) #works
    parser.add_argument("--telo_length", required=True)

    #    parser.add_argument("--restriction_digest_analysis", required=False, default="test")
    return parser

def main(args):
    # convert parameters from strings to usable data types

    args.vrr_length = args.vrr_length == "true"
    args.telo_length = args.telo_length == "true"
    vrr_dict = {}
    telo_length = {}

    out_df = []

    # open stats files
    telo_fh = open("sample_stats.txt", "w")
    vrr_fh = open("sample_stats.VRR.txt", "w")

    if args.telo_length:
        telo_fh.write("Sample_ID\tNumber_of_Reads\tMean_Telomere_Length\tStandard_Deviation_Telo_Length\tQ1\tQ2\tQ3\tMin_Telo_Length\tMax_Telo_Length\n")

    if args.vrr_length:
        vrr_fh.write("Sample_ID\tNumber_of_Reads\tMean_VRR_Telomere_Length\tStandard_Deviation_VRR_Length\tQ1\tQ2\tQ3\tMin_VRR_Telo_Length\tMax_VRR_Telo_Length\n")

    # for every sample in demultiplexed data write out telomere statistics
    for file in args.stat_files:
        df = pd.read_table(file, sep="\t")
        # stats
        df["sample"] = file.strip().split("/")[-1].split(".")[0]

        if args.telo_length:
            telo_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['telo_length']), \
                                                            np.mean(df['telo_length']), np.std(df['telo_length']), np.quantile(df['telo_length'], 0.25), np.quantile(df['telo_length'], 0.5), \
                                                                np.quantile(df['telo_length'], 0.75), min(df['telo_length']), max(df['telo_length'])))

        if args.vrr_length:
            vrr_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['vrr_telo_length']), \
                                                                np.mean(df['vrr_telo_length']), np.std(df['vrr_telo_length']), np.quantile(df['vrr_telo_length'], 0.25), np.quantile(df['vrr_telo_length'], 0.5), \
                                                                    np.quantile(df['vrr_telo_length'], 0.75), min(df['vrr_telo_length']), max(df['vrr_telo_length'])))

    telo_fh.close()
    vrr_fh.close()

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
