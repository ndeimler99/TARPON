#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat_files", nargs="+", required=True)
    parser.add_argument("--vrr_length", required=True) #works

    #    parser.add_argument("--restriction_digest_analysis", required=False, default="test")
    return parser

def main(args):

    vrr_dict = {}
    telo_length = {}

    out_fh = open("sample_stats.txt", "w")
    out_df = []

    if args.vrr_length:
        vrr_fh = open("sample_stats.VRR.txt", "w")
        vrr_fh.write("Sample_ID\tNumber_of_Reads\tMean VRR Telomere Length\tQ1\tQ2\tQ3\tMin VRR Telo Length\tMax VRR Telo Length\n")

    out_fh.write("Sample_ID\tNumber_of_Reads\tMean_Telomere_Length\tStandard Deviation Telo Length\tQ1\tQ2\tQ3\tMin Telo Length\tMax Telo Length\n")

    for file in args.stat_files:
        print(file)
        df = pd.read_table(file, sep="\t")
        # stats
        df["sample"] = file.strip().split("/")[-1].split(".")[0]
        out_df.append(df)

        out_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['telo_length']), \
                                                            np.mean(df['telo_length']), np.std(df['telo_length']), np.quantile(df['telo_length'], 0.25), np.quantile(df['telo_length'], 0.5), \
                                                                np.quantile(df['telo_length'], 0.75), min(df['telo_length']), max(df['telo_length'])))
        
        telo_length[file.strip().split("/")[-1].split(".")[0]] = list(df['telo_length'])
        
        if args.vrr_length:
            vrr_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['vrr_telo_length']), \
                                                                np.mean(df['vrr_telo_length']), np.std(df['vrr_telo_length']), np.quantile(df['vrr_telo_length'], 0.25), np.quantile(df['vrr_telo_length'], 0.5), \
                                                                    np.quantile(df['vrr_telo_length'], 0.75), min(df['vrr_telo_length']), max(df['vrr_telo_length'])))

            vrr_dict[file.strip().split("/")[-1].split(".")[0]] = list(df['vrr_telo_length'])
            
    combined_pd = pd.concat(out_df)
    combined_pd.to_csv("combined_df.csv", index=False)

    out_fh.close()
    if args.vrr_length:
        vrr_fh.close()

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)