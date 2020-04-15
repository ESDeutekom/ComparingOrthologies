#python3
import sys
from itertools import combinations, combinations_with_replacement

if len(sys.argv) != 3:
    print("Need 2 arguments: [Parsed BioGRID interactions file] [pseudo negative interaciton out file name]")
    sys.exit()

biogrid_file = sys.argv[1] # file parsed with all non redundant intractions biogrid
pseudo_negatives = sys.argv[2] #out file with all pseudo negatives

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()

#From the non-redundat interactions parsed from BioGrid
#Retrieve all the sequences that have no interactions with eachother
#But are found in the set at least X times

biogrid_file = open(biogrid_file, "r")

int_set = {} #interactions from biogrid
#make biogrid interactions dictionary
for lines in biogrid_file:
    line = lines.rstrip().split()
    Euk4_A = line[2]
    Euk4_B = line[5]
    if len(Euk4_A.split(",")) >1 or len(Euk4_B.split(",")) >1: # only take the ones where there are not multiple IDs
        pass
    else:
        if Euk4_A not in int_set:
            int_set[Euk4_A] = {}
            int_set[Euk4_A][Euk4_B] = [line]
        else:
            if Euk4_B not in int_set[Euk4_A]:
                int_set[Euk4_A][Euk4_B] = [line]
            else:
                print("somthing is wrong", Euk4_A, Euk4_B)
biogrid_file.close()

biogrid_file = sys.argv[1]
biogrid_file = open(biogrid_file, "r")
#Counts of interactions for individual interactors
#We can use this later to select which ones to keep as non interacting pair
Interactors_d = {}
for lines in biogrid_file:
    line = lines.rstrip().split()
    Euk4_A = line[2]
    Euk4_B = line[5]
    if len(Euk4_A.split(",")) >1 or len(Euk4_B.split(",")) >1:
        pass
    else:
        if Euk4_A not in Interactors_d:
            Interactors_d[Euk4_A] = 1
        else:
            Interactors_d[Euk4_A] += 1
        if Euk4_B not in Interactors_d:
            Interactors_d[Euk4_B] = 1
        else:
            Interactors_d[Euk4_B] += 1

#select only interactors that have more than X times been in an interaction
x = 5
new_Interactors = {}
for key, value in Interactors_d.items():
    if value > x:
        new_Interactors[key] = value

pseudo_negatives = open(pseudo_negatives, "w")
#get all possible pairs from these interactors
all_pairs_dict = {k: 1 for k in combinations_with_replacement(list(new_Interactors.keys()), 2)}
#see which ones are (not) interacting
for pair in all_pairs_dict:
    #if pair[0] == pair[1]:
        #print(pair)
    if pair[0] in int_set:
        if pair[1] not in int_set[pair[0]]: #so not interacting
            pseudo_negatives.write("\t".join(["\t".join(pair), "\n"]))
        else:
            print(pair)
    elif pair[1] in int_set:
        if pair[0] not in int_set[pair[1]]:
            pseudo_negatives.write("\t".join(["\t".join(pair), "\n"]))
        else:
            print(pair)
