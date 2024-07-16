#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

input_files = [sys.argv[1]]
vrr_length = sys.argv[2] == "true"

vrr_dict = {}
telo_length = {}

out_fh = open("sample_stats.txt", "w")
out_df = []

if vrr_length:
    vrr_fh = open("sample_stats.VRR.txt", "w")
    vrr_fh.write("Sample_ID\tNumber of Reads\tMean VRR Telomere Length\tQ1\tQ2\tQ3\tMin VRR Telo Length\tMax VRR Telo Length\n")

out_fh.write("Sample_ID\tNumber of Reads\tMean Telomere Length\tStandard Deviation Telo Length\tQ1\tQ2\tQ3\tMin Telo Length\tMax Telo Length\n")

for file in input_files:
    df = pd.read_table(file, sep="\t")
    # stats
    df["sample"] = file.strip().split("/")[-1].split(".")[0]
    out_df.append(df)

    for i in df['telo_length']:
        if i == 0:
            print(i)
    out_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['telo_length']), \
                                                           np.mean(df['telo_length']), np.std(df['telo_length']), np.quantile(df['telo_length'], 0.25), np.quantile(df['telo_length'], 0.5), \
                                                            np.quantile(df['telo_length'], 0.75), min(df['telo_length']), max(df['telo_length'])))
    
    telo_length[file.strip().split("/")[-1].split(".")[0]] = list(df['telo_length'])
    
    if vrr_length:
        vrr_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file.strip().split("/")[-1].split(".")[0], len(df['vrr_telo_length']), \
                                                               np.mean(df['vrr_telo_length']), np.std(df['vrr_telo_length']), np.quantile(df['vrr_telo_length'], 0.25), np.quantile(df['vrr_telo_length'], 0.5), \
                                                                np.quantile(df['vrr_telo_length'], 0.75), min(df['vrr_telo_length']), max(df['vrr_telo_length'])))

        vrr_dict[file.strip().split("/")[-1].split(".")[0]] = list(df['vrr_telo_length'])
        
combined_pd = pd.concat(out_df)
combined_pd.to_csv("combined_df.csv", index=False)

out_fh.close()
if vrr_length:
    vrr_fh.close()