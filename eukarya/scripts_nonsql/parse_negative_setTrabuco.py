#python3
###################################
## Parse negative interaction set
## Trabuco et al. 2012
###################################
import sys
from collections import Counter


if len(sys.argv) != 3:
    print("Need 2 arguments: [Trabuco negative interaction file input] [uniprot to enseml id] [metadata human proteome] [out file name]")
    sys.exit()

interaction_file = sys.argv[1] #negative interaction file of Trabuco et al. 2012
ukb_ens_file = sys.argv[2] #human uniprotkb id to ensembl id
euk4_file = sys.argv[3] #metadata file of our human proteome
negative_file = sys.argv[4] #out file

try:
	open(sys.argv[1])
    open(sys.argv[2])
    open(sys.argv[3])
except IOError:
    print("No such input file"); sys.exit()

#translate uniprot kb to ensembl
ukb_to_ens = {}
for lines in open(ukb_ens_file):
    line = lines.rstrip().split("\t")
    ens = line[1]
    if ens == "Ensembl" or ens == "EnsemblGenome": #human or yeast respectivly
        ukb_id = line[0]
        ens_id = line[2]
        if ukb_id in ukb_to_ens:
            ukb_to_ens[ukb_id] += [ens_id]
        else:
            ukb_to_ens[ukb_id] = [ens_id]

euk4_file = open(euk4_file, "r")
euk4_file.readline() #skip first line
#Only one protein per ensemble ID
#translate enseml to eukarya
ens_to_euk = {}
for lines in euk4_file:
    line = lines.rstrip().split("\t")
    seqID = line[0]
    L_t = line[7]
    #ENSG = line[9][:line[9].find(".")] #gives -1 if not found, this is dangerous of you don't have a dot, because then line[9][:-1] is taken
    #This fortunatly does not matter if is is really ensembl IDs, because they always have dot in our data for hsap
    #for scer this is not the case
    ENSG = line[9].split(".")[0]
    if L_t == "1":
        ens_to_euk[ENSG] = seqID
euk4_file.close()


# Look at the negative set and translate to eukarya 4
nInt_d = {}
for lines in open(interaction_file, "r"):
    line = lines.rstrip().split("\t")
    nIntA = line[0].split(":")[1]
    nIntB = line[1].split(":")[1]
    shortest_p = str(line[14].split(":")[1])
    pair = (nIntA, nIntB)
    if pair[0] in ukb_to_ens:
        if pair[1] in ukb_to_ens:
            ens1 = ukb_to_ens[pair[0]] #can be list of ensemble IDs
            ens2 =  ukb_to_ens[pair[1]]
            for ens1e in ens1: #for every ensemble id in the list
                for ens2e in ens2:
                    if ens1e in ens_to_euk:
                        #print(ens1e)
                        if ens2e in ens_to_euk:
                            if pair in nInt_d:
                                #translate uniprotkb into ensembl
                                # Every ensemble should only have one eukarya4
                                #so every ensembl pair should give one eukarya4 pair
                                nInt_d[pair][(ens1e, ens2e)] = (ens_to_euk[ens1e], ens_to_euk[ens2e])
                            else:
                                nInt_d[pair] = {}
                                nInt_d[pair][(ens1e, ens2e)] = (ens_to_euk[ens1e], ens_to_euk[ens2e])

#The same unprot pairs can have different ensemble and hsap euk4 pairs
header = ["UnipkbA", "UnipkbB", "ENS_A", "ENS_B", "Euk4_A", "Euk4_B","\n"]
negative_file = open(negative_file, "w")
negative_file.write("\t".join(header))

for key, value in nInt_d.items():
    for key2, value2, in value.items():
        negative_file.write("\t".join(["\t".join(key), "\t".join(key2), "\t".join(value2), "\n"]))
