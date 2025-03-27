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
    parser.add_argument("--repeat_distribution", required=True)
    
    return parser

def hamming_distance(seq, repeat):
    hamming = 0
    for i in range(0, len(seq)):
        if seq[i] != repeat[i]:
            hamming += 1
    return hamming

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
        

def main(args):

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
    wt_count = 0
    telo_sequences = []
    with open(args.stats_out, "w") as stats_out_fh:
        # create new stats file with percentage of WT repeat, % of mutant repeat, % of other one nucleotide variations within each read - I can than histogram this in R and html report - as well as create bar plot and sort bar plot by vrr length
        stats_out_fh.write("read_id\tstrand\tread_len\tvrr_start_pos\tvrr_telo_length\ttelo_length\tread_qual\ttelo_qual\twt_composition\tone_nucl_variant_composition\n")
        for aln in aln_file:
            print(aln.query_name)
            telo_seq = aln.query_sequence[int(telo_dict[aln.query_name][3]):]
            
            wt_nucl = telo_seq.count(args.repeat) * len(args.repeat)
            variant_nucl = len(list(regex.finditer(r"(%s){s<=1}" % args.repeat, telo_seq))) * len(args.repeat)

            wt_comp = wt_nucl / len(telo_seq) * 100
            variant_comp = variant_nucl / len(telo_seq) * 100 - wt_comp
            telo_dict[aln.query_name].extend([str(wt_comp), str(variant_comp)])
            stats_out_fh.write("\t".join(telo_dict[aln.query_name]) + "\n")
            telo_sequences.append(telo_seq)

            wt_count += telo_seq.count(args.repeat)
            telo_seq = telo_seq.replace(args.repeat, "N"*len(args.repeat))
            for rep in one_nucl_dict:
                one_nucl_dict[rep] += telo_seq.count(rep)
                telo_seq = telo_seq.replace(rep, "N"*len(rep))
            print(wt_count)
    
    print(one_nucl_dict)
    for rep in one_nucl_dict:
        one_nucl_dict[rep] = one_nucl_dict[rep]/(sum(one_nucl_dict.values()) + wt_count) * 100
    
    with open(args.repeat_distribution, "w") as repeat_fh:
        for rep in one_nucl_dict:
            for i in range(0, len(args.repeat)):
                if rep[i] != args.repeat[i]:
                    repeat_fh.write("{}\t{}\t{}\t{}\n".format(rep, i, rep[i], one_nucl_dict[rep]))

    with open("telo_sequences.txt", "w") as out_fh:
        for seq in telo_sequences:
            out_fh.write(seq + "\n")

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)