#!/hosts/linuxhome/scarab/eva2/Programs/miniconda3/bin/python
#python3
import sys
from itertools import combinations
from euk_og_dict import *

from scipy.spatial import distance

################################################
#get inter orthology distances of profiles
################################################

#Run from main directory
og_files = ["./eukarya/annotations/Orthogroups_broccoli.txt",\
"./eukarya/annotations/Orthogroups_eggnog_diamond.txt",\
"./eukarya/annotations/Orthogroups_eggnog_hmmer_corrected.txt",\
"./eukarya/annotations/Orthogroups_orthofinder_blast_e-3.txt",\
"./eukarya/annotations/Orthogroups_orthofinder_diamond_e-3.txt",\
"./eukarya/annotations/Orthogroups_panther_different.txt",\
"./eukarya/annotations/Orthogroups_Sonicparanoid_sensitive.txt",\
"./eukarya/annotations/Orthogroups_Swiftortho_c50.txt"]

#Because the orthofinder blast and swift ortho has two species missing, all profiles for the orthers must be with the same 2 removed
#Other wise you get different sized profiles when comparing orthofinder to the rest
og_profiles = ["./Results/Profiles/broccoli_binary",\
"./Results/Profiles/eggnog_diamond_binary",\
"./Results/Profiles/eggnog_hmmer_corrected_binary",\
"./Results/Profiles/orthofinder_blast_e-3_binary",\
"./Results/Profiles/orthofinder_diamond_e-3_binary",\
"./Results/Profiles/panther_different_binary",\
"./Results/Profiles/Sonicparanoid_sensitive_binary",\
"./Results/Profiles/Swiftortho_c50_binary"]

leca_files = ["./eukarya/annotations/leca_orthologous_group_list_broccoli",\
"./eukarya/annotations/leca_orthologous_group_list_eggnog_diamond",\
"./eukarya/annotations/leca_orthologous_group_list_eggnog_hmmer_corrected",\
"./eukarya/annotations/leca_orthologous_group_list_orthofinder_blast_e-3",\
"./eukarya/annotations/leca_orthologous_group_list_orthofinder_diamond_e-3",\
"./eukarya/annotations/leca_orthologous_group_list_panther_different",\
"./eukarya/annotations/leca_orthologous_group_list_Sonicparanoid_sensitive",\
"./eukarya/annotations/leca_orthologous_group_list_Swiftortho_c50"]

out_dir = './Results/All-vs-all-profiles/'
#get index combinations to calculate all-vs-all inter get_distances
combis = list(combinations(range(len(og_files)),2))

#get euk to og translations for orthologies and their profiles
for (i,j) in combis:
    name1 = leca_files[i].split("/")[-1].split("leca_orthologous_group_list_")[1]
    name2 = leca_files[j].split("/")[-1].split("leca_orthologous_group_list_")[1]
    new_name = name1 + "_vs_" + name2
    print(new_name)
    euk_to_og1, profile_dict1 = euk_to_og(og_files[i], og_profiles[i], leca_files[i])
    euk_to_og2, profile_dict2 = euk_to_og(og_files[j], og_profiles[j], leca_files[j])
    #compare the ortholgies and translate them to each other and get profiles
    og_dist_d = {}
    for euk in euk_to_og1:
        if euk in euk_to_og2:
            #multiple sequences could have the same OGs --> Giving the same distances,
            #check if they are not already in distance dict og_dist_d
            og1 = euk_to_og1[euk]
            og2 = euk_to_og2[euk]
            if (og1, og2) not in og_dist_d or (og2, og1) not in og_dist_d:
                p1 = profile_dict1[og1]
                p2 = profile_dict2[og2]
                #if one of the profiles is larger, remove EGRA/NUSP from profile
                #index are 55 and 93
                if len(p1) > len(p2): #different amount of species
                    p1.pop(93) #NUSP
                    p1.pop(55) #EGRA
                elif len(p2) > len(p1): #different amount of species
                    p2.pop(93) #NUSP
                    p2.pop(55) #EGRA
                og_dist_d[(og1, og2)] = distance.cosine(p1, p2)
    out_file = "".join([out_dir,new_name])
    out_file = open(out_file, "w")
    for key, value in og_dist_d.items():
        out_file.write("\t".join([",".join(key), str(round(value, 2)), "\n"]))
