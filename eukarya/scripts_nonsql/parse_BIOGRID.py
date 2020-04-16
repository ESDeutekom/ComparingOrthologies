#python3

import sys
from itertools import groupby, chain
from operator import itemgetter, add
import copy

#################################
## Parse BIOGRID interactions
#################################

if len(sys.argv) != 8:
    print("Need 7 arguments: [BioGRID interactions input file] [metadatafile of set] [ID translation file] [all_out file name] [sys_out file name] [taxID] [minimum pubmed IDs wanted]")
    sys.exit()

BIOGRID_file = sys.argv[1] #interactions from BIOGRID database
Euk4_file = sys.argv[2] #metadatafile of taxID
IDtranslation_file = sys.argv[3] #translation file between ensembl id's used in biogrid and our set
all_out_interactions = sys.argv[4] #interactions out file for human on our set with
sys_out_interactions = sys.argv[5] #same as above, only also check amount of experimental methods used
taxID = sys.argv[6] #taxID
parse_num=sys.argv[7] #number of pubID's minimal

try:
	open(sys.argv[1])
    open(sys.argv[2])
    open(sys.argv[3])
except IOError:
    print("No such input file(s)"); sys.exit()

#Example taxID
#HSAP_tax = 9606
#DMEL_tax = 7227
#SCER_tax = 559292

#Eukarya file, translation from ensembl to sequencese
def ens_to_euk(euk_meta):
    euk_m = open(euk_meta, "r")
    next(euk_m) #skip first line
    euk_to_ens = {}
    for lines in euk_m:
        line = lines.rstrip().split("\t")
        seqID = line[0]
        #ENSG = line[9][:line[9].find(".")] #gives -1 if not found, this is dangerous of you don't have a dot, because then line[9][:-1] is taken
        #This fortunatly does not matter if is is really ensembl IDs, because they always have dot in our data for hsap
        #for scer this is not the case
        ENSG = line[9].split(".")[0]
        longest_trans = line[7]
        if longest_trans == "1": #only take the longest transcript
            if ENSG in euk_to_ens:
                if seqID not in euk_to_ens[ENSG]:
                    euk_to_ens[ENSG] += [seqID]
            else:
                euk_to_ens[ENSG] = [seqID]
    return euk_to_ens

#ID translation dictionary, not all entrez have ensembl
def enz_to_ens(tax_ID, translation_file):
    translation_file = open(translation_file, "r")
    header = translation_file.readline()
    IDtranslation_dict = {}
    for lines in translation_file:
        line = lines.rstrip().split()
        taxID = line[0] #Organism ID NCBI Taxonomy ID
        if tax_ID == taxID:
            if len(line) > 2:
                Entz = line[1] #GeneID == EtrezGene
                ENSG = line[2] # Ensembl ID
                if str(taxID) == str(tax_ID):
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
    BIOGRID_file = open(BIOGRID_file, "r")
    header = BIOGRID_file.readline()
    for lines in BIOGRID_file:
        line = lines.rstrip().split("\t")
        line_count +=1
        #BIOGRID interaction ID. Is a bit useless, because of redundant interactions having different IDs
        BG_int_ID = line[0]
        #Entrez IDs for interactors A and B & BIOGRID IDs for interactors A and B
        Entz_A = line[1]
        BG_A = line[3]
        Entz_B = line[2]
        BG_B = line[4]
        pubID = line[14] #Publication ID
        taxID = line[15] #Organism ID NCBI Taxonomy ID
        #Info we want to know to check how much we can trust the data
        Sys = line[11] #System
        SysT = line[12] #Sytem type
        TP = line[17] #Throughput
        if str(taxID) == str(tax_ID): # check if it is the right taxID
            if Entz_A in inter_dict:
                if Entz_B in inter_dict[Entz_A]:
                    if pubID in inter_dict[Entz_A][Entz_B]:
                        for i in range(len(inter_dict[Entz_A][Entz_B][pubID])):
                            if Sys not in inter_dict[Entz_A][Entz_B][pubID][i]:
                                inter_dict[Entz_A][Entz_B][pubID] += [Sys]
                    else:
                        inter_dict[Entz_A][Entz_B][pubID] = [Sys]
                else:
                    inter_dict[Entz_A][Entz_B] = {}
                    inter_dict[Entz_A][Entz_B][pubID] = [Sys]
            else:
                    inter_dict[Entz_A] = {}
                    inter_dict[Entz_A][Entz_B] = {}
                    inter_dict[Entz_A][Entz_B][pubID] = [Sys]
    print("lines in BioGRID file: ", line_count)
    return inter_dict

#clean all the empty parts in the nested dictionary
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
                        for sys in BAd[pubID]:
                            if sys in ABd[pubID]:
                                index_sys = nonR_dict[ABr[0]][ABr[1]][pubID].index(sys)
                                nonR_dict[ABr[0]][ABr[1]][pubID][index_sys] = ""
                            else:
                                nonR_dict[AB[0]][AB[1]][pubID].append(sys)
                                index_sys = nonR_dict[ABr[0]][ABr[1]][pubID].index(sys)
                                nonR_dict[ABr[0]][ABr[1]][pubID][index_sys] = ""
                    else:
                        nonR_dict[AB[0]][AB[1]][pubID] = BAd[pubID]
                        nonR_dict[ABr[0]][ABr[1]][pubID] = ""
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
nonR_dict = nonR_dict(interactionD) #[EntzA][EntzB][pubID] #all non redundant interactions
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

interactions_out_sys = open(sys_out_interactions, "w")
interactions_out_sys.write(header + "\n")
interactions_out_sys.close()

count_int = 0 #count interactions with more than 2 publications
count_int_sys = 0#count interactions with more than 2 experimental validations and publications
print("Going into nonR_dict")
pubID_num = int(parse_num)

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
                        ",".join(list(nonR_dict[A][B].keys()))])
                        if len(list(nonR_dict[A][B].keys())) >= pubID_num: #more than pubID_num publications
                            interactions_out = open ( all_out_interactions,"r" )
                            lineList = interactions_out.readlines()[-1] #readlast line
                            interactions_out.close()
                            if lineList != new_line:
                                interactions_out = open(all_out_interactions,"a")
                                interactions_out.write(new_line+"\n")
                                interactions_out.close()
                                count_int+= 1
                            for i in range(len(list(nonR_dict[A][B].values()))): #if more than pubID_num publications
                                if len(list(nonR_dict[A][B].values())[i]) > 1: #atleast need 2 different experimental validations
                                    interactions_out_sys = open ( sys_out_interactions,"r" )
                                    lineList_sys = interactions_out_sys.readlines()[-1].rstrip()
                                    interactions_out_sys.close()
                                    if lineList_sys != new_line:
                                        interactions_out_sys = open(sys_out_interactions,"a")
                                        interactions_out_sys.write(new_line+"\n")
                                        interactions_out_sys.close()
                                        count_int_sys += 1

print("Interactions translated to eukarya with more than", pubID_num, "publication: ", count_int)
print("Interactions translated to eukarya with more than 1 experimental validation: ", count_int_sys)
