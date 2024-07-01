import gzip
import sys

input_file = sys.argv[1]
repeat = sys.argv[2]
repeat_count = int(sys.argv[3])
c_strand_only = sys.argv[4] == 'true'
out_file = sys.argv[5]
non_telo = sys.argv[6]


def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])


if c_strand_only:
    print("It Worked")
    print(c_strand_only)
    print(type(c_strand_only))

with gzip.open(input_file, 'rt') as input_file_fh, open(out_file, 'w') as out_fh, open(non_telo, "w") as non_telo_fh:
    linecount = 0
    read = []
    for line in input_file_fh:
        linecount += 1
        read.append(line.strip())
        if linecount % 4 == 0:
            if c_strand_only and read[1].count(rev_complement(repeat)) >= repeat_count:
                out_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3])) 
            elif read[1].count(repeat) > repeat_count or read[1].count(rev_complement(repeat)) >= repeat_count:
                out_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3]))    
            else:
                non_telo_fh.write('{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3]))
            read = []
