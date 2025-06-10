#!/usr/bin/env python3

import gzip
import argparse
import pysam

def main(args):

    args.min_percentage = float(args.min_percentage)
    header = ""
    stats_dict = {}
    with open(args.stats_fh, "r") as stats:
        linecount = 1
        for line in stats:
            if linecount == 1:
                linecount += 1
                header = line.strip()
                continue
            line = line.strip().split()
            stats_dict[line[0]] = line
    
    cluster_dict = {}
    cluster_sizes = {}
    with open(args.cluster_results, 'r') as cluster_fh:
        linecount = 1
        for line in cluster_fh:
            if linecount == 1:
                linecount += 1
                continue
            line = line.strip().split("\t")
            count = 0
            for read_id in line[-1].split(","):
                count += 1
                cluster_dict[read_id] = line[3]
            cluster_sizes[line[3]] = count
    
    removed_clusters = []
    for cluster in cluster_sizes:
        if cluster_sizes[cluster] <= sum(cluster_sizes.values()) * args.min_percentage/100:
            removed_clusters.append(cluster)
    
    for read in stats_dict:
        if read in cluster_dict:
            if cluster_dict[read] in removed_clusters:
                cluster_dict[read] = "NA"
        else:
            cluster_dict[read] = "NA"
        stats_dict[read].append(cluster_dict[read])

    with open(args.new_stats_fh, "w") as stats_fh:
        stats_fh.write(header + "\tCluster\n")
        for read in stats_dict:
            stats_fh.write("\t".join(stats_dict[read]))
            stats_fh.write("\n")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_fh", required=True)
    parser.add_argument("--cluster_results", required=True)
    parser.add_argument("--new_stats_fh", required=True)
    parser.add_argument("--min_percentage", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
