#!/hosts/linuxhome/scarab/eva2/Programs/miniconda3/bin/python
#python3
import sys

#Translate ortholgies to find inter orthology differences with distances later on
#Calculate distances between profiles of the same sequences/ogs from different orthologies
#using distance that came up the best in our initial analysis

#Get HSAP sequences from the different orthologies and make translation dict
def euk_to_og(og_file, profile_file, leca_file):
    leca_file = open(leca_file, "r")
    leca_d = {}
    #make a leca dictionary for easy access
    for lines in leca_file:
        leca_id = lines.rstrip()
        leca_d[leca_id] = "-"
    #dictionary for translating euk sequence to og sequences
    euk_og_dict = {}
    og_file = open(og_file, "r")
    for lines in og_file:
        line = lines.rstrip().split(":")
        og_id = line[0]
        seqs = line[1].split(" ")
        HSAP_seqs = [seq for seq in seqs if "HSAP" in seq]
        for HSAP in HSAP_seqs:
            if og_id in leca_d:
                euk_og_dict[HSAP] = og_id
    #make a profile dictionary
    profile_file = open(profile_file, "r")
    profile_file.readline().split("\t")#skip_header
    profile_dict = {}
    for lines in profile_file:
        line = lines.rstrip().split("\t")
        og_id = line[0]
        profile = line[1:len(line)]
        profile = [int(num) for num in profile]
        if og_id in leca_d:
            profile_dict[og_id] = profile
    return euk_og_dict, profile_dict
