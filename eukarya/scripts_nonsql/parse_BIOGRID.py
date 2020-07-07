#python3

import os
import sys
from itertools import groupby, chain
from operator import itemgetter, add
import copy

#######################################################
## Parse BIOGRID interactions for biogrid redundancy
## Translate to our eukarya4 species set
#######################################################

if len(sys.argv) != 6:
    print("Need 5 arguments: [BioGRID interactions input file] [metadatafile of set] [ID translation file] [all_out file name] [taxID]")
    sys.exit()

BIOGRID_file = sys.argv[1] #interactions from BIOGRID database
Euk4_file = sys.argv[2] #metadatafile of taxID
IDtranslation_file = sys.argv[3] #translation file between ensembl id's used in biogrid and our set
all_out_interactions = sys.argv[4] #interactions out file for human on our set with
taxID = sys.argv[5] #taxID

try:
    open(sys.argv[1])
    open(sys.argv[2])
    open(sys.argv[3])
except IOError:
    print("No such input file(s)"); sys.exit()

#Check if file is not empty
for file in (sys.argv[1], sys.argv[2], sys.argv[3]):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

#Example taxID
#HSAP_tax = 9606
#DMEL_tax = 7227
#SCER_tax = 559292

#Eukarya file, translation from ensembl to sequencese
def ens_to_euk(euk_meta):
    euk_m = open(euk_meta, "r")
    next(euk_m) #skip first line
    ens_to_euk = {}
    for lines in euk_m:
        line = lines.rstrip().split("\t")
        seqID = line[0]
        ENSG = line[9].split(".")[0]
        longest_trans = line[7]
        if longest_trans == "1": #only take the longest transcript
            if ENSG in ens_to_euk:
                if seqID not in ens_to_euk[ENSG]:
                    ens_to_euk[ENSG] += [seqID]
            else:
                ens_to_euk[ENSG] = [seqID]
    return ens_to_euk

#ID translation dictionary, not all entrez have ensembl
#some ensembl have more than 1 entrez
def enz_to_ens(tax_ID, translation_file):
    translation_file = open(translation_file, "r")
    header = translation_file.readline()
    IDtranslation_dict = {}
    for lines in translation_file:
        line = lines.rstrip().split()
        taxID = line[0] #Organism ID NCBI Taxonomy ID
        if tax_ID == taxID:
            if len(line) > 2: #if there is a translation
                Entz = line[1] #GeneID == EtrezGene
                ENSG = line[2] # Ensembl ID

                if Entz in IDtranslation_dict:
                    if ENSG not in IDtranslation_dict[Entz]:
                        IDtranslation_dict[Entz] += [ENSG]
                else:
                    IDtranslation_dict[Entz] = [ENSG]
    return IDtranslation_dict

#Cristeria: must be non-redundanrt interactions
#BIOGRID file contains interactions A-->B and B-->A.
#However, they can differ in type of system/experiment etc. and publications
#Since we want to know all publications, we need to check for this
#And not just say the interaction is redundant if they are bidirectional.

#Interaction dictionary
def interaction_dict(tax_ID, BIOGRID_file):
    inter_dict = {}
    line_count = 0
    BIOGRID_file = open(BIOGRID_file, "r") #BioGrid file is tab3 format
    header = BIOGRID_file.readline()
    for lines in BIOGRID_file:
        line = lines.rstrip().split("\t")
        line_count +=1
        #BIOGRID interaction ID. Is a bit useless, because of redundant (same) interactions having different IDs
        BG_int_ID = line[0]
        #Entrez IDs for interactors A and B & BIOGRID IDs for interactors A and B
        Entz_A = line[1]
        BG_A = line[3]
        Entz_B = line[2]
        BG_B = line[4]
        pubID = line[14] #Publication ID
        taxID = line[15] #Organism ID NCBI Taxonomy ID
        #Info we want to know to check how much we can trust the data
        if str(taxID) == str(tax_ID): # check if it is the right taxID
            if Entz_A in inter_dict:
                if Entz_B in inter_dict[Entz_A]:
                    if pubID not in inter_dict[Entz_A][Entz_B]:
                        inter_dict[Entz_A][Entz_B] += [pubID]
                else:
                    inter_dict[Entz_A][Entz_B] = [pubID]
            else:
                inter_dict[Entz_A] = {}
                inter_dict[Entz_A][Entz_B] = [pubID]


    print("lines in BioGRID file: ", line_count)
    return inter_dict

#function to clean all the empty parts in the nested dictionary
def clean_empty(d):
    if not isinstance(d, (dict, list)):
        return d
    if isinstance(d, list):
        return [v for v in (clean_empty(v) for v in d) if v]
    return {k: v for k, v in ((k, clean_empty(v)) for k, v in d.items()) if v}

#then see if the reverse of the pairs B--> A is later in the list
def nonR_dict(interaction_dict):
    #Note:This removes the directional data from BIOGRID
    #remove all redundant interactions from dictionary by tracking
    # A-->B and B--> A with same publication
    #First save all the pairs A-->B
    nonR_dict = copy.deepcopy(interaction_dict) #dictionary to adjust, deep copy of original
    pair_list = [] #interaction pairs
    for A, value in interaction_dict.items():
        for B in value:
            pair_list.append([A,B])
    #make pairs and their reverse
    for i in range(len(pair_list)):
        AB = pair_list[i] #unique pairs but bidirectional
        ABr = AB[::-1] #the reverse of AB
        if  ABr in pair_list: #see if the reverse of AB is in the list
            ir = pair_list.index(ABr) #get the index of the reverse
            #Because reverse of A-->B is B-->A, but reverse of B-->A is again A-->B
            #going through the list, i must always be smaller than ir
            if i < ir: #we need to check if they are the same or what difference is
                #Check if it is the same pubID
                ABd = interaction_dict[AB[0]][AB[1]]
                BAd = interaction_dict[ABr[0]][ABr[1]]
                for pubID in BAd:
                    if pubID in ABd:
                        index_pub = nonR_dict[ABr[0]][ABr[1]].index(pubID)
                        nonR_dict[ABr[0]][ABr[1]][index_pub] = "" #remove
                #else if pubID not in A-B, we want to remove B-A and add pubID from B-A to A-B
                #this way we only keep the A-B non redundant ones
                    else:
                        nonR_dict[AB[0]][AB[1]].append(pubID) #add
                        index_pub = nonR_dict[ABr[0]][ABr[1]].index(pubID)
                        nonR_dict[ABr[0]][ABr[1]][index_pub] = "" #remove
    return(clean_empty(nonR_dict))

# get all non redundant interactions with more than x pubIDs
# get all non redundant interactions with more than sys in pubID as resque?

#translation dictionaries
enz_to_ens = enz_to_ens(taxID, IDtranslation_file) #[Entz] = [ENSG, ENSG, ...]
ens_to_euk = ens_to_euk(Euk4_file) #[ENSG] = [euk4, euk4, ...]

print("Genes in metadata with Ensemble ID: ", len(ens_to_euk.keys()))
print("Entrez from biogrid translated to ensemble: ", len(enz_to_ens.keys()))

#initialize
print("starting")
interactionD = interaction_dict(taxID, BIOGRID_file) #contains all interactions, also redundant ones
print("Interaction BioGRID done")
nonR_dict = nonR_dict(interactionD) #[EntzA][EntzB]: [pubID, pubID] #all non redundant interactions
print("Non redundant interactions done")


R_interactions = 0
for A, Bd in interactionD.items():
    R_interactions += len(list(Bd.keys()))
print("Interactions biogrid: ", R_interactions)
nR_interactions = 0
for A, Bd in nonR_dict.items():
    nR_interactions += len(list(Bd.keys()))
print("Non redundant interactions biogrid: ", nR_interactions)

#write non redundant interactions to file
header = "\t".join(["Entz_A", "ENS_A", "Euk4_A", "Entz_B", "ENS_B", "Euk4_B", "pubIDs"])

interactions_out = open(all_out_interactions,"w")
interactions_out.write(header + "\n")
interactions_out.close()


count_int = 0 #count interactions with more than 2 publications
count_int_sys = 0#count interactions with more than 2 experimental validations and publications
print("Going into nonR_dict")

for A in nonR_dict: # interactor A entrez from biogrid
    for B in nonR_dict[A]: # interactor B entrez from biogrid
        if A in enz_to_ens: # is entrez A translatable to ensemble
            if B in enz_to_ens: # is entrez A translatable to ensemble
                starterA = [] #could be more ensebles to eukarya id's (not a lot, but some)
                starterB = []
                for ensA in enz_to_ens[A]: #then translate to ensemble
                    for ensB in enz_to_ens[B]:
                        ########print(ensA, ensB)
                        if ensA in ens_to_euk: #is ensemble A in eukarya 4
                            if ensB in ens_to_euk:  #is ensemble B in eukarya 4
                                if ensA not in starterA:  #add the to list
                                    starterA +=  ens_to_euk[ensA]
                                if ensB not in starterB:
                                    starterB += ens_to_euk[ensB]
                if starterA: #if the list is not empty for A
                    if starterB: #if the list is not empty for B
                        #then we make a newline to write to OG_file
                        new_line = "\t".join([
                        A, ",".join(enz_to_ens[A]), ",".join(starterA),
                        B, ",".join(enz_to_ens[B]), ",".join(starterB),
                        ",".join(nonR_dict[A][B])])
                        #if len(list(nonR_dict[A][B].keys())) >= pubID_num: #more than pubID_num publications
                        interactions_out = open ( all_out_interactions,"r" )
                        lineList = interactions_out.readlines()[-1] #readlast line
                        interactions_out.close()
                        if lineList != new_line:
                            interactions_out = open(all_out_interactions,"a")
                            interactions_out.write(new_line+"\n")
                            interactions_out.close()
                            count_int+= 1


print("Interactions translated to eukarya: ", count_int)
