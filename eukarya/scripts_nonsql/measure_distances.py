#python3
############################################################
## Measure distances between phylogenetic profiles
## between negative sets and the positive interaction set
############################################################

import sys
import random
import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.stats import spearmanr, pearsonr, kendalltau
from itertools import combinations, combinations_with_replacement

if len(sys.argv) != 6:
    print("Need 5 arguments: [Ortholgoy profiles] [LECA list file] [BioGRID interactions] [Trabuco negative interactions] [Pseudo negative interactions]")
    sys.exit()

pro_file = sys.argv[1] #phylogenetic profiles
leca_file = sys.argv[2] #leca obtained from dollo
int_file = sys.argv[3] #interactions file
RnegSet_file = sys.argv[4] #Trabuco set
PnegSet_file =  sys.argv[5] #pseudo negative set
out_file = sys.argv[6] #file containing all the distances between sets

try:
	open(sys.argv[1])
    open(sys.argv[2])
    open(sys.argv[3])
    open(sys.argv[4])
    open(sys.argv[5])
except IOError:
    print("No such input file"); sys.exit()


#leca file to parse out only the LECA OGs
leca_file = open(leca_file, "r")
leca_og_d = {}
for lines in leca_file:
    leca = lines.rstrip()
    leca_og_d[leca] = 0
leca_file.close()
print("Read leca file length: ", len(leca_og_d) )

def get_distances(pro_file, leca_og_d ,int_file, out_file, funcs, type):
    pro_file = open(pro_file, "r") #OG profile file from different ortholgies
    header= pro_file.readline()
    for lines in pro_file:
        line = lines.rstrip().split("\t")
        og_id = line[0]
        if og_id in leca_og_d:
            profile = [float(i) for i in line[1:]]
            leca_og_d[og_id] = profile #Get profiles for LECA OGs
    pro_file.close()
    print("Read OG file lenght: ", len(leca_og_d))

    int_file = open(int_file, "r") #Interactions or non interactions from biogrid or negative sets
    int_dict = {}
    for lines in int_file:
        line = lines.rstrip().split("\t")
        OG1 = line[-2]
        OG2 = line[-1]
        int_dict[(OG1, OG2)] = True
    int_file.close()
    print("Read interactions length: ", len(int_dict))

    #get all possible combinations of pair of OGs. All possible combinations of OGs in leca og dictionary
    pairs_list = combinations_with_replacement(list(leca_og_d.keys()), 2)

    i_distance_matrix = {} #interacting distance matrix

    for pair in pairs_list: #get the profiles of each OG
        OG1 = pair[0]
        OG2 = pair[1]
        profA = leca_og_d[OG1]
        profB = leca_og_d[OG2]
        for measure, func in funcs.items():
            if (OG1, OG2) in int_dict or (OG2, OG1) in int_dict: #since the niteractions could be from A-->B or B-->A or combinations_with_replacement
            #we need to check both, but don't want redundancy
                if pair in i_distance_matrix:
                    if measure not in i_distance_matrix[pair]:
                        if measure in ["spearman", "pearson","kendalltau"]:
                            dist_measure = func(profA, profB)
                            i_distance_matrix[pair][measure] =  1-dist_measure[0]
                        else:
                            dist_measure = func(profA, profB)
                            i_distance_matrix[pair][measure] =  round(dist_measure,4)
                else:
                    i_distance_matrix[pair] = {}
                    if measure in ["spearman", "pearson", "kendalltau"]:
                        dist_measure = func(profA, profB)
                        i_distance_matrix[pair][measure] =  1-dist_measure[0]
                    else:
                        dist_measure = func(profA, profB)
                        i_distance_matrix[pair][measure] =  round(dist_measure,4)

    print("Got interacting matrix length: ", len(i_distance_matrix))#, len(ni_distance_matrix))
    out_file = open(out_file, "a")

    for row in i_distance_matrix:
        new_line = "\t".join([",".join(row), "\t".join([str(i_distance_matrix[row][i]) for i in sorted(i_distance_matrix[row])]), type, "\n"])
        out_file.write(new_line)
    out_file.close()

#Dice gives all 0's for some reason
#Correlation gives NaN values --> Give profiles entering correltation a noise
#distances to calculate
funcs = {"euclidean": distance.euclidean, "braycurtis": distance.braycurtis, \
"cityblock": distance.cityblock, "cosine": distance.cosine, "jaccard": distance.jaccard, \
"dice": distance.dice, "yule": distance.yule, "kulsinski": distance.kulsinski,
"russellrao": distance.russellrao, "sokalmichener": distance.sokalmichener,
"dice": distance.dice, "rogerstanimoto": distance.rogerstanimoto, "spearman": spearmanr,"kendalltau": kendalltau}


header = ['pair'] + sorted(list(funcs.keys()))+ ["Interaction"]

out_files = open(out_file, "w")
out_files.write("\t".join(header)+"\n")
out_files.close()

print("Biogrid distances")
get_distances(pro_file, leca_og_d, int_file, out_file, funcs, type = "BioGrid")
print("RusselNeg distances")
get_distances(pro_file, leca_og_d, RnegSet_file, out_file, funcs, type = "RusselNeg")
print("PseudoNeg distances")
get_distances(pro_file, leca_og_d, PnegSet_file, out_file, funcs, type = "PseudoNeg")
