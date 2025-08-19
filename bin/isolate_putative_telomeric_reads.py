#!/usr/bin/env python3

import gzip
import argparse
import pysam
import multiprocessing

def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def isolate_reads(sequence_dict, repeat, repeat_count, c_strand_only, mutant):

    telo_reads = []
    non_telo_reads = []

    for read in sequence_dict:
        if mutant == "false":
            if c_strand_only and sequence_dict[read].count(rev_complement(repeat)) >= repeat_count:
                telo_reads.append(read)
            elif not c_strand_only and (sequence_dict[read].count(repeat) >= repeat_count or sequence_dict[read].count(rev_complement(repeat)) >= repeat_count):
                telo_reads.append(read)
            else:
                non_telo_reads.append(read)
        else:
            if c_strand_only and sequence_dict[read].count(rev_complement(repeat)) + sequence_dict[read].count(rev_complement(mutant)) >= repeat_count:
                telo_reads.append(read)
            elif not c_strand_only and (sequence_dict[read].count(repeat) + sequence_dict[read].count(mutant) >= repeat_count or sequence_dict[read].count(rev_complement(repeat)) + sequence_dict[read].count(rev_complement(mutant)) >= repeat_count):
                telo_reads.append(read)
            else:
                non_telo_reads.append(read)
    return telo_reads, non_telo_reads

def split_dict(d, n):
    """Split dictionary d into n smaller dictionaries."""
    items = list(d.items())
    k, m = divmod(len(items), n)
    chunks = []
    for i in range(n):
        start = i*k + min(i, m)
        end = (i+1)*k + min(i+1, m)
        chunks.append(dict(items[start:end]))
    return chunks

def main(args):

    start = time.time()
    # convert parameters from strings to usable data types

    args.c_strand_only = args.c_strand_only == "true"
    args.repeat_count = int(args.repeat_count)
    args.threads = int(args.threads)

    # for read in gzipped fastq file perform analysis
    # if the argument c strand only is set convert the telomeric repeat is C strand and identify reads
    # otherwise check for both forward and reverse telomeric repeats
    input_file_fh = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    out_fh = pysam.AlignmentFile(args.out_file, "wb", template=input_file_fh)
    non_telo_fh = pysam.AlignmentFile(args.non_telo, "wb", template=input_file_fh)

    sequence_dict = {}
    aln_dict = {}
    for aln in input_file_fh:
        aln.query_name = aln.query_name.split()[0]
        sequence_dict[aln.query_name] = aln.query_sequence
        aln_dict[aln.query_name] = aln

    bam_time = time.time()
    bam_chunks = split_dict(sequence_dict, args.threads)  # split into 4 equal sized chunks

    mp_args = [(chunk, args.repeat, args.repeat_count, args.c_strand_only, args.mutant) for chunk in bam_chunks]
    
    # for i, c in enumerate(chunks):
    #     print(f"Chunk {i} has {len(c)} reads")
    with multiprocessing.Pool(args.threads) as pool:
        results = pool.starmap(isolate_reads, mp_args)
    
    telo_reads = [read for out in results for read in out[0]]
    non_telo_reads = [read for out in results for read in out[1]]
    
    for read in telo_reads:
        out_fh.write(aln_dict[read])
    
    for read in non_telo_reads:
        non_telo_fh.write(aln_dict[read])
    input_file_fh.close()
    out_fh.close()
    non_telo_fh.close()

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",required=True)
    parser.add_argument("--repeat",required=True)
    parser.add_argument("--repeat_count", required=True)
    parser.add_argument("--c_strand_only", required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--non_telo", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--threads", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)