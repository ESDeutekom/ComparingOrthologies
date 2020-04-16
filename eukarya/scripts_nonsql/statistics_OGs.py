#python3

import sys
import statistics as s

##########################
## statistics on orthogroups
##########################

if len(sys.argv) != 3:
    print("Need 2 arguments: [Orthologous group input file] [LECA orthologous groups input list]")
    sys.exit()

OG_file = sys.argv[1] #Orthoglous groups file
leca_file = sys.argv[2] #LECA orthologous groups files

try:
	open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()


leca_file = open(leca_file, "r")
leca_og_d = {}
for line in leca_file:
    og_id = line.rstrip()
    leca_og_d[og_id] = True

def counts_OG(OG_file, leca_og_d):
    count_dict = {} #contains all statistics
    og_dict = {} #to count per og species
    count_dict["count_OGs_all"] = 0 # count the amount of OGs inferred from method
    count_dict["count_singlets"] = 0 #all OGs with only 1 protein
    count_dict["count_groups"] = 0 #count OG with more than 1 protein
    seq_counts = [] #list for mean and median with sequence counts per OG
    seq_leca_counts = [] #sequence counts for OGs in LECA
    count_species_per_OG = [] #unique species per OG counting
    total_prots_assigned = 0 #total proteins assigned to an OG, but a real group and not a singlet
    total_prots = 0 #total proteins assigned in the dataset
    leca_og = 0 #counts of LECA OGs
    OG_open = open(OG_file, "r")
    for lines in OG_open:
        line = lines.rstrip().split(":")
        OG_id = line[0]
        orgsL = line[1].split()
        total_prots += len(orgsL)

        count_dict["count_OGs_all"] += 1 #every line is a OG from method (singlets and rest)
        if OG_id in leca_og_d:
            seq_leca_counts += [len(orgsL)]
            leca_og += 1
        if len(orgsL) > 1: # a real OG
            og_dict[OG_id] = []
            count_dict["count_groups"] += 1
            seq_counts += [len(orgsL)]
            total_prots_assigned += len(orgsL)
            for org in orgsL:
                org_id = org[0:4]
                if org_id not in og_dict[OG_id]:
                    og_dict[OG_id] += [org_id] #count per OG the # of species

        else: #single protein OGs, not real "group"
            count_dict["count_singlets"] +=1
    for key, values in og_dict.items():
        count_species_per_OG += [len(values)]
    max_species = max(count_species_per_OG)
    count_dict["median_og_size"] = s.median(seq_counts)
    count_dict["mean_og_size"] = round(s.mean(seq_counts),1)
    count_dict["max_og_size"] = max(seq_counts)
    count_dict["min_og_size"] = min(seq_counts)
    count_dict["single_species_OGs"] = count_species_per_OG.count(1)
    count_dict[" ".join(["OGs with all", str(max_species), "species present"])] = count_species_per_OG.count(max_species)
    count_dict["Total proteins assigned"] = total_prots_assigned
    count_dict["Total proteins"] = total_prots
    count_dict["Percentage assigned proteins"] = round((float(total_prots_assigned)/float(total_prots))*100, 1)
    count_dict["LECA OGs"] = leca_og
    count_dict["median_leca_size"] = s.median(seq_leca_counts)
    count_dict["mean_leca_size"] = round(s.mean(seq_leca_counts),1)
    count_dict["stdev_leca_size"] = round(s.stdev(seq_leca_counts),1)
    count_dict["max_leca_size"] = max(seq_leca_counts)
    count_dict["min_leca_size"] = min(seq_leca_counts)

    return count_dict

dict_out = counts_OG(OG_file, leca_og_d)
for key, value in dict_out.items():
    print (key, value) #print to file
