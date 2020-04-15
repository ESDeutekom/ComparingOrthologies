#python3
###############################################
## Get a BioGrid corrected version of
## Trabuco et al. negative interaction neg_set
###############################################

import sys

if len(sys.argv) != 3:
    print("Need 2 arguments: [BioGRID interactions parsed] [negative interaction Trabuco et al. parsed]")
    sys.exit()

biogrid_file = sys.argv[1] #BioGRID interactions parsed
neg_set_file = sys.argv[2] #Trabuco et al. negative set parsed
overlap_out = sys.argv[3] #interactions that are overlapping
neg_set_crosscheck = sys.argv[4] #interaction that are not overlapping (and this negative)

try:
	open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()


biogrid_file = open(biogrid_file, "r")
header = biogrid_file.readline()
print(header)
pairs_list_bio = {}
for lines in biogrid_file:
    line = lines.rstrip().split()
    Euk4_A = line[2]
    Euk4_B = line[5]
    pairs_list_bio[(Euk4_A, Euk4_B)] = "True"
    pairs_list_bio[(Euk4_B, Euk4_A)] = "True"
biogrid_file.close()

neg_set_file = open(neg_set_file, "r")
overlap_out = open(overlap_out, "w")
neg_set_crosscheck = open(neg_set_crosscheck, "w")
header2 = neg_set_file.readline()
print(header2)
neg_set_crosscheck.write(header2)
for lines in neg_set_file:
    line = lines.rstrip().split()
    nEuk4_A = line[4]
    nEuk4_B = line[5]
    if (nEuk4_A, nEuk4_B) in  pairs_list_bio:
        overlap_out.write(lines)
    elif  (nEuk4_B, nEuk4_A) in  pairs_list_bio:
        overlap_out.write(lines)
    else:
        neg_set_crosscheck.write(lines)
neg_set_file.close()
overlap_out.close()
neg_set_crosscheck.close()
