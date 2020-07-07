#python 3

import os
import sys

######################################################
## Make required OG input file from orthology output
## Broccoli orthology
######################################################
try:
    open(sys.argv[1])
except IOError:
    print("No such input file(s)"); sys.exit()

#Check if file is not empty
if os.path.getsize(sys.argv[1]) <= 1:
    print(sys.argv[1], "file is empty"); sys.exit()

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
