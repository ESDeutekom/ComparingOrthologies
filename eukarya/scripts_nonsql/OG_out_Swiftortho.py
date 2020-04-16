#python 3

import sys

######################################################
## Make required OG input file from orthology output
## SiwftOrtho orthology
######################################################

file_og = sys.argv[1]
file_sequence_ids = sys.argv[2]
file_out = sys.argv[3]

sequence_id_dict = {}
for line in open(file_sequence_ids, 'r'):
    line = line.rstrip().split(": ")
    id = line[0].split("_")
    id_ = "|".join(id)
    sequence = line[1]
    sequence_id_dict[id_] = sequence

line_count = 1
file_out = open(file_out, 'w')

for line in open(file_og, 'r'):
    new_line = []
    seq_line = line.rstrip().split("\t")
    for seqID in seq_line:
        new_line += [sequence_id_dict[seqID]]
    new_string = " ".join(new_line)
    file_out.write("%s: %s\n" % (str(line_count), new_string))
    line_count += 1
