#python3

import os
import sys
import random
import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.stats import spearmanr, pearsonr, kendalltau
from itertools import combinations, combinations_with_replacement

#######################################################
## Measure distances for the random interaction set
## This was done later, so added seperate file and code
#######################################################
if len(sys.argv) != 4:
    print("Need 3 arguments: [Orthologous profile input file] [random interactions file] [distance out file name]")
    sys.exit()

pro_file = sys.argv[1] #phylogenetic profile
int_file = sys.argv[2] #random interactions file
out_file = sys.argv[3] #out file with distances between random set

try:
    open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

def get_distances(pro_file, int_file, out_file, funcs, type):
    pro_file = open(pro_file, "r") #OG profile file from different ortholgies
    header= pro_file.readline()

    og_d = {} #make dictionary of OGs and their profiles
    for lines in pro_file:
        line = lines.rstrip().split("\t")
        og_id = line[0]
        profile = [float(i) for i in line[1:]]
        og_d[og_id] = profile #Get profiles for LECA OGs
    pro_file.close()

    int_file = open(int_file, "r") #get the pairs from the random pairs file
    random_dict = {}
    for lines in int_file:
        line = lines.rstrip().split(",")
        OG1 = line[0]
        OG2 = line[1]
        random_dict[(OG1, OG2)] = True
    int_file.close()

    i_distance_matrix = {} #interacting distance matrix
    for pair in random_dict.keys(): #get the profiles of each OG
        OG1 = pair[0]
        OG2 = pair[1]
        profA = og_d[OG1]
        profB = og_d[OG2]
        for measure, func in funcs.items():
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

# \"sokalsneath": distance.sokalsneath,  "pearson": pearsonr,
header = ['pair'] + sorted(list(funcs.keys()))+ ["Interaction"]

out_files = open(out_file, "w")
out_files.write("\t".join(header)+"\n")
out_files.close()

get_distances(pro_file, int_file, out_file, funcs, type = "RandomSet")
