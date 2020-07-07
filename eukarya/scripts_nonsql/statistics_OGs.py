#python3

import os
import sys
import statistics as s
import pandas as pd
from find_loss import *
##########################
## statistics on orthogroups
##########################

"""if len(sys.argv) != 5:
print("Need 4 arguments: [Orthologous group input file] [LECA orthologous groups input list] [method name] [stats out file name]")
sys.exit()
"""

OG_file = sys.argv[1] #Orthoglous groups file
leca_file = sys.argv[2] #LECA orthologous groups files
met_name = sys.argv[3] #name of orthology
dollo_tree = sys.argv[4]
out_file = sys.argv[5] #out file

try:
    open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

for file in (sys.argv[1], sys.argv[1]):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

leca_file1 = open(leca_file, "r")
leca_og_d = {}
for line in leca_file1:
    og_id = line.rstrip()
    leca_og_d[og_id] = True

print(len(leca_og_d))

def counts_OG(OG_file, leca_og_d):

    total_165=2811230
    total_167=2865661

    count_dict = {} #contains all statistics
    og_dict = {} #to count per og species
    #count_dict["Number of annotated OGs (incl. singlets)"] = 0 # count the amount of OGs inferred from method
    #count_dict["count single protein OGs"] = 0 #all OGs with only 1 protein
    count_dict["Number of OGs"] = 0 #count OG with more than 1 protein

    seq_counts = [] #list for mean and median with sequence counts per OG
    seq_leca_counts = [] #sequence counts for OGs in LECA

    count_species_per_OG = [] #unique species per OG counting

    total_seqs_assigned = 0 #total proteins assigned to an OG, but a real group and not a singlet
    total_seqs = 0 #total proteins assigned in the dataset

    leca_og = 0 #counts of LECA OGs
    OG_open = open(OG_file, "r")

    for lines in OG_open:
        line = lines.rstrip().split(":")
        OG_id = line[0]
        orgsL = line[1].split()
        total_seqs += len(orgsL) #total proteins in orthology annotated to OG (singlet or not)

        #count_dict["Number of annotated OGs (incl. singlets)"] += 1 #every line is a OG from method (singlets and rest)
        if OG_id in leca_og_d:
            seq_leca_counts += [len(orgsL)] #total proteins in orthology annoted to LECA OG
            leca_og += 1 #total leca_ogs

        if len(orgsL) > 1: # a real OG (aka group of more than 1 sequence)
            og_dict[OG_id] = []
            count_dict["Number of OGs"] += 1
            seq_counts += [len(orgsL)] #counts of sequences in og (for median and mean)
            total_seqs_assigned += len(orgsL)
            for org in orgsL:
                org_id = org[0:4]
                if org_id not in og_dict[OG_id]:
                    og_dict[OG_id] += [org_id] #count per OG the # of species

        #else: #single protein OGs, not real "group"
        #    count_dict["count single protein OGs"] +=1
    for key, values in og_dict.items():
        count_species_per_OG += [len(values)]
    max_species = max(count_species_per_OG)
    count_dict["Median OG size"] = s.median(seq_counts)
    count_dict["Mean OG size"] = round(s.mean(seq_counts),1)
    #count_dict["max OG size"] = max(seq_counts)
    #count_dict["min OG size"] = min(seq_counts)
    #count_dict["single species OGs"] = count_species_per_OG.count(1)
    #count_dict[" ".join(["OGs with all", str(max_species), "species present"])] = count_species_per_OG.count(max_species)
    if max_species == 165:
        count_dict["% proteins assigned by orthology"] = round((float(total_seqs)/float(total_165))*100,1)
        count_dict["% proteins assigned to LECA OG from total"] = round((float(sum(seq_leca_counts))/float(total_165))*100,1)
    if max_species == 167:
        count_dict["% proteins assigned by orthology"] = round((float(total_seqs)/float(total_167))*100,1)
        count_dict["% proteins assigned to LECA OG from total"] = round((float(sum(seq_leca_counts))/float(total_167))*100,1)
    #count_dict["Total proteins"] = total_seqs
    count_dict["% assigned proteins to OGs"] = round((float(total_seqs_assigned)/float(total_seqs))*100, 1)

    count_dict["Number LECA OGs"] = leca_og
    count_dict["Median LECA OG size"] = s.median(seq_leca_counts)
    count_dict["Mean LECA OG size"] = round(s.mean(seq_leca_counts),1)
    count_dict["stdev LECA OG size"] = round(s.stdev(seq_leca_counts),1)
    count_dict["Max LECA OG size"] = max(seq_leca_counts)
    count_dict["% to LECA OG assigned proteins"] =round((float(sum(seq_leca_counts))/float(total_seqs_assigned))*100,1)
    return count_dict

dict_out = counts_OG(OG_file, leca_og_d)

loss_dict,_ = loss_dict(dollo_tree, leca_file) #returns loss_dict and independent loss distributions, only need loss dict
dict_out.update(loss_dict)
dict_df = pd.DataFrame.from_dict(dict_out, orient='index', columns = [str(met_name)])
print(dict_df)
if os.path.exists(out_file):
    df = pd.read_csv(out_file, sep = ",", index_col = 0)
    df_out = pd.concat([df, dict_df], axis=1)#add this to already existing files/calculations
    df_out.to_csv(out_file)#, header = False)"""
else:
    dict_df.to_csv(out_file)
