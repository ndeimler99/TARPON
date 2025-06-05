#!/usr/bin/env python3

import regex
import argparse
import pysam

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--stats_out", required=True)
    parser.add_argument("--wt_processivity", required=True)
    parser.add_argument("--mt_processivity", required=True)
    parser.add_argument("--mutant", required=True)

    return parser


def calculate_variant_nucl(telo_seq, repeat, mutant):
    
    variant_nucl = len(list(regex.finditer(r"(%s){s<=1}" % repeat, telo_seq))) * len(args.repeat)
    if hamming_distance(repeat, mutant) <= 1:
        variant_nucl = variant_nucl - telo_seq.count(mutant) * len(mutant)
    
    return variant_nucl

def update_processivity_dict(processivity_dict, telo_seq, repeat, mutant):

    wt_split = telo_seq.split(repeat)
    mt_split = telo_seq.split(mutant)

    for subseq in wt_split[1:(len(wt_split)-1)]:
        if subseq == "":
            continue
        mt_count = subseq.count(mutant)
        if mt_count * len(mutant) >= 0.8 * len(subseq):
            if mt_count in processivity_dict["mutant"]:
                processivity_dict["mutant"][mt_count] += 1
            else:
                processivity_dict["mutant"][mt_count] = 1
    
    for subseq in mt_split[1:(len(mt_split)-1)]:
        if subseq == "":
            continue
        wt_count = subseq.count(repeat)
        if wt_count * len(repeat) >= 0.8*len(subseq):
            if wt_count in processivity_dict["wild_type"]:
                processivity_dict["wild_type"][wt_count] += 1
            else:
                processivity_dict["wild_type"][wt_count] = 1
    return processivity_dict

def main(args):
    telo_sequences = []
    telo_dict = {}
    #line[3] = telo_start
    #line[1] = strand
    with open(args.stats_file, "r") as stats_fh:
        linecount = 0 
        for line in stats_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            telo_dict[line[0]] = line

    aln_file = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    # loop through reads

    wt_count = 0
    mt_count = 0

    processivity_dict = {"wild_type":{}, "mutant":{}}

    with open(args.stats_out, "w") as stats_out_fh:
        # create new stats file with percentage of WT repeat, % of mutant repeat, % of other one nucleotide variations within each read - I can than histogram this in R and html report - as well as create bar plot and sort bar plot by vrr length
        
        stats_out_fh.write("read_id\ttelo_length\twt_perc\tmt_perc\tother_perc\n")

        for aln in aln_file:
            telo_seq = aln.query_sequence[int(telo_dict[aln.query_name][3]):]
            telo_sequences.append(telo_seq)
            wt_nucl = telo_seq.count(args.repeat) * len(args.repeat)
            mt_nucl = telo_seq.count(args.mutant) * len(args.mutant)

            mt_comp = mt_nucl / len(telo_seq) * 100
            wt_comp = wt_nucl / len(telo_seq) * 100

            telo_dict[aln.query_name].extend([str(wt_comp), str(mt_comp)])
            stats_out_fh.write("{}\t{}\t{}\t{}\n".format(aln.query_name, telo_dict[aln.query_name][4], wt_comp, mt_comp, 100-wt_comp-mt_comp))

            # calculate processivity
            processivity_dict = update_processivity_dict(processivity_dict, telo_seq[-1000:], args.repeat, args.mutant)

    with open("telo_sequences.txt", "w") as out_fh:
        for seq in telo_sequences:
            out_fh.write(seq + "\n")

    with open(args.wt_processivity, "w") as wt_processivity_fh:
        for i in range(1, max(processivity_dict["wild_type"])):
            if i in processivity_dict["wild_type"]:
                wt_processivity_fh.write("{}\t{}\n".format(i, processivity_dict["wild_type"][i]))
            else:
                wt_processivity_fh.write("{}\t0\n".format(i))

    with open(args.mt_processivity, "w") as mt_processivity_fh:
        for i in range(1, max(processivity_dict["mutant"])):
            if i in processivity_dict["mutant"]:
                mt_processivity_fh.write("{}\t{}\n".format(i, processivity_dict["mutant"][i]))
            else:
                mt_processivity_fh.write("{}\t0\n".format(i))

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)