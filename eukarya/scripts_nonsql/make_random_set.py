#python3

import os
import sys
import random
####################################################################
## Make a random set of OGs (random pairs) to compare with the
## interacting and non-interacting pairs and their distances
## we need to get a random subset
####################################################################

if len(sys.argv) != 4:
    print("Need 3 arguments: [Orthologous group input file] [LECA list file] [random set out file name]")
    sys.exit()

OG_file = sys.argv[1]
LECA_file = sys.argv[2]
out_file = sys.argv[3]

try:
    open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

#Check if file is not empty
for file in (sys.argv[1], sys.argv[2]):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

LECA_dict = {}
LECA_file = open(LECA_file, "r")
for line in LECA_file:
    LECA = line.rstrip()
    LECA_dict[LECA] = True
print("leca size: ", len(LECA_dict))
LECA_file.close()

OG_list = []
OG_file = open(OG_file, 'r')

for line in OG_file:
    line = line.rstrip().split(': ')
    OG = line[0]
    org_list = line[1].split(' ')
    if any('HSAP' in sequence for sequence in org_list): #OG must contain human
        if OG in LECA_dict: #OG must be in LECA (otherwise you get distances between animal vs plant only etc.)
            OG_list += [OG]
OG_file.close()

#To make it most fair to compare to the other sets:
#also do not take OG pairs that are the same (no phylogenetic value)
#do not have redundancy (A--> B, but not B--> A)

random_dict={}
print(len(OG_list)," is amount of OGs containing HSAP")

while len(random_dict) < len(OG_list):
    OG_pair = (random.sample(OG_list, 1)[0], random.sample(OG_list, 1)[0]) #from the set, randomly select an OG
    if OG_pair not in random_dict:
        if OG_pair[::-1] not in random_dict:
            if OG_pair[0] != OG_pair[1]:
                random_dict[OG_pair] = True

out_file = open(out_file,"w")
for key in random_dict:
    out_file.write(",".join([key[0], key[1]])+'\n')
