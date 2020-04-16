#python3

import sys
from random import sample

####################################################################
## Make a random set of OGs (random pairs) to compare with the
## interacting and non-interacting pairs and their distances
## we need to get a random subset
####################################################################

if len(sys.argv) != 4:
    print("Need 2 arguments: [Orthologous group input file] [LECA list file] [random set out file name]")
    sys.exit()

OG_file = sys.argv[1]
LECA_file = sys.argv[2]
out_file = sys.argv[3]

try:
	open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

LECA_dict = {}
LECA_file = open(LECA_file, "r")
for line in LECA_file:
    LECA = line.rstrip()
    LECA_dict[LECA] = True
print(len(LECA_dict))
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

A_random = sample(OG_list, len(OG_list)) #from the set, randomly select an OG
B_random = sample(OG_list, len(OG_list))

print(len(OG_list))
OG_pairs = zip(A_random, B_random) #add them together to form a (random) pair
out_file = open(out_file,"w")
for el in list(OG_pairs):
    out_file.write(",".join([el[0], el[1]])+'\n')
