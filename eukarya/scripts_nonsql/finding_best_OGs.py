#python3

import os
import sys
from itertools import combinations, product
import random
import pandas as pd
import numpy as np
############################################################
## Find overlap between (the best) OGs
## Rev. comment: Can we find which OGs are very comparable between all the methods?
############################################################
out_file = sys.argv[1]
out_file2 = sys.argv[2]

OG_FILES=["./eukarya/annotations/Orthogroups_orthofinder_blast_e-3.txt",\
"./eukarya/annotations/Orthogroups_orthofinder_diamond_e-3.txt",\
"./eukarya/annotations/Orthogroups_broccoli.txt",\
"./eukarya/annotations/Orthogroups_eggnog_diamond.txt",\
"./eukarya/annotations/Orthogroups_eggnog_hmmer_corrected.txt",\
"./eukarya/annotations/Orthogroups_panther_different.txt",\
"./eukarya/annotations/Orthogroups_Sonicparanoid_sensitive.txt",\
"./eukarya/annotations/Orthogroups_Swiftortho_c50.txt"]

LECA_FILES=["./eukarya/annotations/leca_orthologous_group_list_orthofinder_blast_e-3",\
"./eukarya/annotations/leca_orthologous_group_list_orthofinder_diamond_e-3",\
"./eukarya/annotations/leca_orthologous_group_list_broccoli",\
"./eukarya/annotations/leca_orthologous_group_list_eggnog_diamond",\
"./eukarya/annotations/leca_orthologous_group_list_eggnog_hmmer_corrected",\
"./eukarya/annotations/leca_orthologous_group_list_panther_different",\
"./eukarya/annotations/leca_orthologous_group_list_Sonicparanoid_sensitive",\
"./eukarya/annotations/leca_orthologous_group_list_Swiftortho_c50"]

#Check if file is not empty
for file in (OG_FILES):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

#Check if file is not empty
for file in (LECA_FILES):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

outfile=sys.argv[1]

def make_leca_og_dict(leca_file, og_file):

    leca_dict = {}
    leca_file = open(leca_file, "r")
    for line in leca_file:
        leca = line.rstrip()
        leca_dict[leca] = True
    leca_file.close()

    og_dict_leca = {}
    og_file = open(og_file, "r")
    for line in og_file:
        line = line.rstrip().split(":")
        OG_id = line[0] #orthologous group id
        speciesL = line[1].split() #seqeunce id list
        if OG_id in leca_dict:
            og_dict_leca[OG_id] = speciesL
    print("leca ogs", len(og_dict_leca))
    og_file.close()
    return og_dict_leca

overall_dict = {}
for i in range(0,len(OG_FILES)):
    name = OG_FILES[i].split("/")[-1].split("Orthogroups_")[-1].split(".")[0]
    print(name)
    og_dict = make_leca_og_dict(LECA_FILES[i], OG_FILES[i])
    overall_dict[name] = og_dict

all_vs_all = combinations(overall_dict.keys(),2)
#print("combinations of orthologies: ", len(list(all_vs_all)))

df = pd.DataFrame()
df2 = pd.DataFrame()
for combi in all_vs_all:
    print(combi)
    counts = 0
    bi_direct = 0
    perfect = 0
    perfect_pair = []
    bi_pair = []
    orthology1 = overall_dict[combi[0]]
    orthology2 = overall_dict[combi[1]]
    og_combi = list(product(list(orthology1.keys()), list(orthology2.keys())))
    for ogs in og_combi:
        og1 = ogs[0]
        og2 = ogs[1]
        #Unfortunatly, need to remove EGRA and NUSP since we are looking at "perfect" overlap
        seqs1 = set([x for x in orthology1[og1] if x[0:4] not in ["EGRA","NUSP"]])
        seqs2 = set([x for x in orthology2[og2] if x[0:4] not in ["EGRA","NUSP"]])
        intersect = seqs1 & seqs2
        #perfect overlap
        if len(intersect)/len(seqs1) == 1 and len(intersect)/len(seqs2) == 1:
            perfect += 1
            perfect_pair += [ogs]
        #bidirectional overlap
        if len(intersect)/len(seqs1) > 0.95 and len(intersect)/len(seqs2) > 0.95:
            bi_direct += 1
            bi_pair += [ogs]
        #how many orthologous groups are (complelety) contained in the other
        #if len(intersect)/len(seqs1) == 1:
        #    counts += 1

    print("perfect overlap: ", perfect)
    print("more than 95% bidirectional overlap: ", bi_direct)
    df_pair = pd.Series(perfect_pair).rename(combi)
    df_bi = pd.Series(bi_pair).rename(combi)
    df = pd.concat([df, df_pair], axis = 1)
    df2 = pd.concat([df2, df_bi], axis = 1)

df.to_csv(out_file, na_rep='')
df2.to_csv(out_file2, na_rep='')
#for all orthologies, find the best overlapping OGs (top 10 or something)
#(higest intersection) between the methods
#count how much from the OG size is overlapping
#compare to other top 10's
