#python3

###################################################
## Translate negative interactions to OG ids
###################################################

import sys

if len(sys.argv) != 4:
    print("Need 3 arguments: [Negative interaction set] [Ortholgous group file] [negative interaction in OGs out file]")
    sys.exit()

neg_set = sys.argv[1] #negative interaction dataset
OG_file = sys.argv[2] #Orthologous group file
neg_out_file = sys.argv[3] #out fle

try:
	open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()


#eukarya to og depends on orthology used.
#Translate eukarya to OGs
og_file = open(OG_file, "r") #.TXT HAS NO HEADER

euk_to_og_d = {}
for lines in og_file:
    line = lines.rstrip().split(":")
    OG_id = line[0]
    seq_list = line[1].strip().split(" ")
    #if OG_id in leca_d:
    for seq in seq_list:
        if seq not in euk_to_og_d:
            euk_to_og_d[seq] = [OG_id]

print("Length euk to og: ", len(euk_to_og_d))

neg_set  = open(neg_set, "r")
neg_out_file = open(neg_out_file, "w")

for lines in neg_set:
    line = lines.rstrip().split('\t')
    eukA = line[-2]
    eukB = line[-1]
    if eukA in euk_to_og_d and eukB in euk_to_og_d:
        line += euk_to_og_d[eukA]
        line +=  euk_to_og_d[eukB]
        line += ["\n"]
        neg_out_file.write("\t".join(line))
    else:
        line += ["-","-","\n"]
        neg_out_file.write("\t".join(line))
