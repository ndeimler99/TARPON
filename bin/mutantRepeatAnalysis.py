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
    parser.add_argument("--repeat_distribution", required=True)
    parser.add_argument("--mutant", required=True)

    return parser

def generate_one_nucleotide_dict(repeat):
    
    dna_list = ["A", "T", "C", "G"]
    one_nucl_dict = {}
    repeat = list(repeat)
    for i in range(0, len(repeat)):
        for nucl in dna_list:
            if repeat[i] == nucl:
                pass
            else:
                new_repeat = list(repeat)
                new_repeat[i] = nucl
                one_nucl_dict[''.join(new_repeat)] = 0
    
    return one_nucl_dict

def hamming_distance(seq, repeat):
    hamming = 0
    for i in range(0, len(seq)):
        if i > len(repeat):
            return hamming
        if seq[i] != repeat[i]:
            hamming += 1
    return hamming


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
    one_nucl_dict = generate_one_nucleotide_dict(args.repeat)
    if args.mutant not in one_nucl_dict:
        one_nucl_dict[args.mutant] = 0

    wt_count = 0
    mt_count = 0

    processivity_dict = {"wild_type":{}, "mutant":{}}

    with open(args.stats_out, "w") as stats_out_fh:
        # create new stats file with percentage of WT repeat, % of mutant repeat, % of other one nucleotide variations within each read - I can than histogram this in R and html report - as well as create bar plot and sort bar plot by vrr length
        stats_out_fh.write("read_id\tstrand\tread_len\tvrr_start_pos\tvrr_telo_length\ttelo_length\tread_qual\ttelo_qual\twt_composition\tmutant_composition\tone_nucl_variant_composition\n")
        for aln in aln_file:
            telo_seq = aln.query_sequence[int(telo_dict[aln.query_name][3]):]
            telo_sequences.append(telo_seq)
            wt_nucl = telo_seq.count(args.repeat) * len(args.repeat)
            mt_nucl = telo_seq.count(args.mutant) * len(args.mutant)
            
            variant_nucl = calculate_variant_nucl(telo_seq, args.repeat, args.mutant)

            mt_comp = mt_nucl / len(telo_seq) * 100
            wt_comp = wt_nucl / len(telo_seq) * 100

            variant_comp = variant_nucl / len(telo_seq) * 100 - wt_comp

            telo_dict[aln.query_name].extend([str(wt_comp), str(mt_comp), str(variant_comp)])
            stats_out_fh.write("\t".join(telo_dict[aln.query_name]) + "\n")

            # calculate processivity
            processivity_dict = update_processivity_dict(processivity_dict, telo_seq, args.repeat, args.mutant)


            mt_read_count = telo_seq.count(args.mutant)
            wt_count += telo_seq.count(args.repeat)
            telo_seq = telo_seq.replace(args.repeat, "N"*len(args.repeat))
            mt_count += mt_read_count
            telo_seq = telo_seq.replace(args.mutant, "N" * len(args.mutant))
            one_nucl_dict[args.mutant] += mt_read_count
            for rep in one_nucl_dict:
                if rep == args.mutant:
                    continue
                one_nucl_dict[rep] += telo_seq.count(rep)
                telo_seq = telo_seq.replace(rep, "N"*len(rep))
    
    for rep in one_nucl_dict:
        one_nucl_dict[rep] = one_nucl_dict[rep]/(sum(one_nucl_dict.values()) + wt_count) * 100
    
    with open(args.repeat_distribution, "w") as repeat_fh:
        for rep in one_nucl_dict:
            if rep == args.mutant and hamming_distance(args.repeat, args.mutant) > 1:
                repeat_fh.write("{}\t{}\t{}\t{}\n".format(rep, "mutant", "mutant", one_nucl_dict[rep]))
                continue
            for i in range(0, len(args.repeat)):
                if rep[i] != args.repeat[i]:
                    repeat_fh.write("{}\t{}\t{}\t{}\n".format(rep, i, rep[i], one_nucl_dict[rep]))

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