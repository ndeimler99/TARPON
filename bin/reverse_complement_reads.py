import sys

input_file = sys.argv[1]
repeat = sys.argv[2]
threshold = float(sys.argv[3])
c_strand_only = sys.argv[4] == 'true'
out_file = sys.argv[5]
g_strand_out = sys.argv[6]
c_strand_out = sys.argv[7]


def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])


with open(input_file, 'r') as input_file_fh, open(out_file, 'w') as out_fh, open(c_strand_out, 'w') as c_out, open(g_strand_out, 'w') as g_out:
    linecount = 0
    read = []
    for line in input_file_fh:
        linecount += 1
        read.append(line.strip())
        if linecount % 4 == 0:
            c = read[1].count(rev_complement(repeat))
            g = read[1].count(repeat)
            if c_strand_only:
                if c/(c+g) >= threshold:
                    # write out to main file and c_strand file
                    out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                    c_out.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
            else:
                if c/(c+g) < threshold and c/(c+g) > 1-threshold:
                    pass
                elif c/(c+g) >= threshold:
                    #write out to main file and c_strand file
                    out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                    c_out.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'C', rev_complement(read[1]), read[2], ''.join([ i for i in read[3][::-1]])))
                elif g/(c+g) >= threshold:
                    #write out to main file and g_strand file
                    out_fh.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'G', read[1], read[2], read[3]))
                    g_out.write('{}\t{}\n{}\n{}\n{}\n'.format(read[0], 'G', read[1], read[2], read[3]))
                
            read = []
