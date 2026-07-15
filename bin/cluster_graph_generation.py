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
    uniq_keys = set()

    for i, cluster in enumerate(clusters):
        cluster_name = f"Concordance_{i}"
        cluster_dict[cluster_name] = {}

        for node in cluster:
            sample, gene = node.split(".", 1)
            uniq_keys.add(sample)

            # Store the FULL node name for graph lookups
            cluster_dict[cluster_name][sample] = node

    uniq_keys = sorted(uniq_keys)


    with open(args.out, "w") as out_fh:
        header = ["Concordance_Cluster"] + uniq_keys
        out_fh.write("\t".join(header) + "\n")

        for cluster_name, members in cluster_dict.items():
            row = [cluster_name]

            for sample in uniq_keys:
                if sample in members:
                    # Write only the gene name
                    row.append(members[sample].split(".", 1)[1])
                else:
                    row.append("NA")

            out_fh.write("\t".join(row) + "\n")


    for cluster_name, members in cluster_dict.items():

        outfile = f"CONCORDANCE_TABLES/{cluster_name}.aln_stats.txt"
        with open(outfile, "w") as fh:
            fh.write("\t" + "\t".join(uniq_keys) + "\n")

            for sampleA in uniq_keys:
                nodeA = members.get(sampleA)
                row = [sampleA]
                for sampleB in uniq_keys:
                    nodeB = members.get(sampleB)
                    if sampleA == sampleB:
                        row.append("100")
                    elif nodeA is None or nodeB is None:
                        row.append("NA")
                    elif G.has_edge(nodeA, nodeB):
                        row.append(f"{G[nodeA][nodeB]['weight']:.2f}")
                    else:
                        row.append("NA")

                fh.write("\t".join(row) + "\n")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aln_files", required=True, nargs="+")
    parser.add_argument("--out", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)