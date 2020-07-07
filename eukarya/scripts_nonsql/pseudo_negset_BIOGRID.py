#python3

import os
import sys
from itertools import combinations, combinations_with_replacement

###################################################################
## get a pseudo negative set from the interaction (BioGRID) set
##################################################################

if len(sys.argv) != 3:
    print("Need 2 arguments: [Parsed BioGRID interactions file] [pseudo negative interaciton out file name]")
    sys.exit()

biogrid_file = sys.argv[1] # file parsed with all non redundant intractions biogrid with OG id's
pseudo_negatives = sys.argv[2] #out file with all pseudo negatives

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()

#Check if file is not empty
if os.path.getsize(sys.argv[1]) <= 1:
    print(biogrid_file, "file is empty"); sys.exit()

#From the non-redundat interactions parsed from BioGrid
#Retrieve all the sequences that have no interactions with eachother
#But are found in the set at least X times

biogrid_file = open(biogrid_file, "r")

#make biogrid interactions dictionary
int_set = {} #interactions from biogrid
#Counts of interactions for individual interactors
#We can use this later to select which ones to keep as non interacting pair
Interactors_d = {}
for lines in biogrid_file:
    line = lines.rstrip().split()
    Euk4_A = line[2]
    Euk4_B = line[5]
    pair = (Euk4_A, Euk4_B) #non-redundant pairs from biogrid, should only contain A-->B and not B-->A
    if len(Euk4_A.split(",")) >1 or len(Euk4_B.split(",")) >1: # only take the ones where there are not multiple IDs
        pass
    else:
        if pair not in int_set:
            int_set[pair] = True
        elif pair[::-1] in int_set:
            print("reverse interaction", pair[::-1]) #This should not happen, euk id only had one ens id
        else:
            print("pair in interaction", pair) #this could happen because entrez could map to same ensembl (not a lot)
        #count interactors and how many interactions they have
        if Euk4_A not in Interactors_d:
            Interactors_d[Euk4_A] = 1
        else:
            Interactors_d[Euk4_A] += 1
        if Euk4_B not in Interactors_d:
            Interactors_d[Euk4_B] = 1
        else:
            Interactors_d[Euk4_B] += 1
biogrid_file.close()

#select only interactors that have more than X times been in an interaction
x = 5
new_Interactors = {}
for interactor, count in Interactors_d.items():
    if count > x:
        new_Interactors[interactor] = count

pseudo_negatives = open(pseudo_negatives, "w")
#get all possible pairs from these interactors with more than 5 interactions
all_pairs_dict = {k: 1 for k in combinations_with_replacement(list(new_Interactors.keys()), 2)} #function gets (A,A), (A,B), but not (B,A)
#from all interactor pairs with more than 5 interactors each
#see which ones are (not) interacting
count = 0
for pair in all_pairs_dict: #for all possible pairs with interactors with more than 5 interactors
    if pair not in int_set: #not an interacting pair
        if pair[::-1] not in int_set: #if reverse pair not interacting pair
            #if pair[0] != pair[1]: #for the negative set we do not want proteins that are the same (since they will have exactly similar profiles)
            pseudo_negatives.write("\t".join(["\t".join(pair), "\n"]))

pseudo_negatives.close()
