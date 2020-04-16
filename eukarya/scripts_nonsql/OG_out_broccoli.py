#python 3

import sys

######################################################
## Make required OG input file from orthology output
## Broccoli orthology
######################################################

file_in = sys.argv[1]
file_out = sys.argv[2]

file_out = open(file_out, 'w')

for line in open(file_in, "r"):
    line = line.rstrip()
    if line.startswith("#"):
        pass
    else:
        OG_line = line.split("\t")
        OG_id = OG_line[0]
        seqs = OG_line[1]
        file_out.write(OG_id+": " + seqs + "\n")
#file_in.close()
file_out.close()
