#python3

import sys
from Bio import SeqIO
from parse_domtbout_func import *

###################################################
## Parse panther LECA hits, best non overlapping
###################################################

# functions to parse and get best non-overlapping hits and pseudo_genes

if len(sys.argv) != 3:
    print("Need 2 arguments: [domtblout hmmsearch] [hits out file name]")
    sys.exit()

panther_dom = sys.argv[1]
all_out = sys.argv[2]

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()

#gets list of hits
panther_dict = parse_domtbout(panther_dom) #parse the file into dictionary

best_hits = {}
#Best hits without overlap. Overlap was too stringent, was getting weird results
for seq, hits in panther_dict.items():
    #get best hit with best e-value (no overlap checking, just best hits)
    #sort list from best to worst sequence e-value
    hits.sort(key = itemgetter(1), reverse=False)
    #Best panther OG hit is hit[0] after sorting
    best_OG = hits[0]
    best_OG_name = best_OG[0]
    if best_OG_name in best_hits:
        best_hits[best_OG_name] += [seq]
    else:
        best_hits[best_OG_name] = [seq]
all_out = open(all_out, 'w')
for OG, sequences in best_hits.items():
    all_out.write('%s\t%s\n' % (OG + ":", " ".join(sequences)))
all_out.close()
