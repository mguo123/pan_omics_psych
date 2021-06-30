'''
python bed_longrange_formatting.py <bed_file> <output_file>
longrange
The longrange track is a bed format-like file type. Each row contains columns from left to right: chromosome, start position (0-based), and end position (not included), interaction target in this format chr2:333-444,55. As an example, interval “chr1:111-222” interacts with interval “chr2:333-444” on a score of 55, we will use following two lines to represent this interaction:

chr1    111 222  chr2:333-444,55
chr2    333 444  chr1:111-222,55

Important:Be sure to make TWO records for a pair of interacting loci, one record for each locus.


# also need to remove lines like this:
chr1 -1      1       chr2:-1-1,inf


'''
import sys, re

bed_file=sys.argv[1]
output_file = sys.argv[2]
with open(bed_file, 'r') as f:
    with open(output_file, 'w') as g:
        for line in f:
            line_arr = line.strip().split('\t')
            if int(line_arr[1])>0:
                line_arr = line_arr[:4]
                line_arr[0]='chr'+line_arr[0]
                line_arr[3]='chr'+line_arr[3]
                new_line = '\t'.join(line_arr)

                # create record but reversed with same score
                second_node = re.split("[,:-]", line_arr[3])
                new_second_node = line_arr[0] + ':' + line_arr[1] + '-' + line_arr[2] + ',' + second_node[-1]
                reverse_line = '\t'.join(second_node[:-1] + [new_second_node])

                g.write(new_line+'\n')
                g.write(reverse_line+'\n')




