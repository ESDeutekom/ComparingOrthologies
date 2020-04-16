#python 3

import sys

######################################################
## Make required OG input file from orthology output
## SonicParanoid orthology
######################################################

file_in = sys.argv[1]
file_out = sys.argv[2]

file_out = open(file_out, 'w')

#SonicParanoid
for line in open(file_in, "r"):
    if line.startswith("group_id"):
        continue
    line = line.rstrip().split('\t')
    new_line = []
    OG_id = line[0]
    seqs = line[1:]
    while "*" in seqs:
        seqs.remove('*')
    for el in seqs:
        seq = el.split(",")
        new_line += seq

    file_out.write(OG_id+": " + " ".join(new_line) + "\n")
