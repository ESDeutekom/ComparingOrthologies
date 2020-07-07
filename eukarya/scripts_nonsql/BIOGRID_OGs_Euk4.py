#python3

import os
import sys

###################################################
## Translate interactions to OG ids of orthology
## Same file as biogrid, just with OG id or "-"
###################################################

if len(sys.argv) != 4:
    print("Need 3 arguments: [bioGrid interactions parsed] [Orthologous group input file] [orthology interactions output file name]")
    sys.exit()

interactions_BioG = sys.argv[1] #positive interactions BioGRID parsed
OG_file = sys.argv[2] #orthologous groups file
bio_out_file = sys.argv[3] #interactions file out

try:
    open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

#Check if file is not empty
for file in (sys.argv[1], sys.argv[2]):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

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
            euk_to_og_d[seq] = OG_id

print("Length euk to og: ", len(euk_to_og_d))

int_set = open(interactions_BioG, "r")
header2=int_set.readline().rstrip().split("\t")

header2 += ["OG_A","OG_B", "\n"]
bio_out_file = open(bio_out_file, "w")
bio_out_file.write("\t".join(header2))
for lines in int_set:
    line = lines.rstrip().split("\t")
    eukA = line[2]
    eukB = line[5]
    pubIDs = line[-1].split(",")
    if len(pubIDs) > 4: #select only interaction with atleast 5 pubIDs
        if eukA in euk_to_og_d and eukB in euk_to_og_d:
            line += [euk_to_og_d[eukA]]
            line += [euk_to_og_d[eukB]]
            line += ["\n"]
            bio_out_file.write("\t".join(line))
        else:
            if eukA in euk_to_og_d:
                line += [euk_to_og_d[eukA],"-", "\n"]
                bio_out_file.write("\t".join(line))
            elif eukB in euk_to_og_d:
                line += ["-",euk_to_og_d[eukB], "\n"]
                bio_out_file.write("\t".join(line))
            else:
                line += ["-","-", "\n"]
                bio_out_file.write("\t".join(line))
