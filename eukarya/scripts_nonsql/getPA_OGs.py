#python3

import sys

#####################################################
## From orthology inference output
## Get the presence and absence profiles
#####################################################

if len(sys.argv) != 3:
    print("Need 2 arguments: [Orthologous group input file] [PA output file name]")
    sys.exit()

OGs_txt = sys.argv[1]
PA_out = sys.argv[2]

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()


OGs_txt = open(OGs_txt, "r")
presences = {}

for lines in OGs_txt:
    line = lines.rstrip().split(":")
    OG_id = line[0] #orthologous group id
    speciesL = line[1].split() #seqeunce id list
    for i in range(0, len(speciesL)):
        species = speciesL[i][0:4] #get the species name from the sequence id
        if OG_id in presences:
            if species not in presences[OG_id]:
                presences[OG_id] += [species]
        else:
            presences[OG_id] = [species]

#write to temporary file
PA_out = open(PA_out, "w")
for OG_id in presences:
    PA_out.write(OG_id + "\t" + "\t".join(presences[OG_id]) + "\n")
PA_out.close()
