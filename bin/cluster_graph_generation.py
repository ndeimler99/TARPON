#!/usr/bin/env python3

import argparse
import networkx as nx
import igraph as ig
import leidenalg as la

def main(args):
   
    G = nx.Graph()
    for file in args.aln_files:
        best_matches_a = {}
        best_matches_b = {}
        with open(file, "r") as fh:
            linecount = 0
            for line in fh:
                if linecount == 0:
                    linecount += 1
                    continue
                line = line.strip().split()
                if float(line[2]) > 90:
                    if line[0] not in best_matches_a:
                        best_matches_a[line[0]] = [line[1], float(line[2])]
                    else:
                        if float(line[2]) > best_matches_a[line[0]][1]:
                            best_matches_a[line[0]] = [line[1], float(line[2])]

                    if line[1] not in best_matches_b:
                        best_matches_b[line[1]] = [line[0], float(line[2])]
                    else:
                        if float(line[2]) > best_matches_b[line[1]][1]:
                            best_matches_b[line[1]] = [line[0], float(line[2])]
                
        for match in best_matches_a:
            G.add_edge(match, best_matches_a[match][0], weight=best_matches_a[match][1])
        for match in best_matches_b:
            G.add_edge(match, best_matches_b[match][0], weight=best_matches_b[match][1])
           


    #clusters = list(nx.connected_components(G)) 

    # map nodes to integers
    nodes = list(G.nodes())
    idx = {n:i for i,n in enumerate(nodes)}

    edges = [(idx[a], idx[b]) for a,b in G.edges()]
    weights = [G[a][b]['weight'] for a,b in G.edges()]

    g = ig.Graph(edges=edges, directed=False)
    g.es['weight'] = weights

    partition = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        weights='weight'
    )

    clusters = [
        [nodes[i] for i in cluster]
        for cluster in partition
    ]

    cluster_dict = {}
    i = 0
    uniq_keys = []
    for cluster in clusters:
        cluster_dict["Concordance_{}".format(i)] = {}
        for sample in cluster:
            uniq_keys.append(sample.split(".")[0])
            cluster_dict["Concordance_{}".format(i)][sample.split(".")[0]] = sample.split(".")[1]
        i+= 1

    uniq_keys = set(uniq_keys)
    print(uniq_keys)

    with open(args.out, "w") as out_fh:
        out_str = "Concordance_Cluster\t"
        for sample in uniq_keys:
            out_str += sample + "\t"
        out_str.strip()
        out_fh.write(out_str + "\n")
        for cluster in cluster_dict:
            out_str = ""
            out_str += cluster + "\t"
            for sample in uniq_keys:
                if sample in cluster_dict[cluster]:
                    out_str += cluster_dict[cluster][sample] + "\t"
                else:
                    out_str += "NA\t"
            out_str.strip()
            out_fh.write(out_str + "\n")



def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aln_files", required=True, nargs="+")
    parser.add_argument("--out", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)