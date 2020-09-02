#!/hosts/linuxhome/scarab/eva2/Programs/miniconda3/bin/python
#python3
import sys
import random
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, auc
from multiprocessing import Process
from ROC_statistics import *
#Get ROC statistics with permutations and bootstrapping
#There is a large difference in datasize between positive and negative set.
import pandas as pd
#Run from main directory
file_dir = "./Results/Distances/"
#Interaction files mapped to orthologies and their distances between phylogenetic profiles
files = ["eggnog_diamond_distances",\
        "eggnog_hmmer_corrected_distances",\
        "orthofinder_diamond_e-3_distances",\
        "orthofinder_blast_e-3_distances",\
        "broccoli_distances",\
        "panther_different_distances",\
        "Sonicparanoid_sensitive_distances",\
        "Swiftortho_c50_distances"]
out_file = "./Results/roc_bootstraps"

file_list = ["".join([file_dir, el]) for el in files]

def bootstrap_all(file_list, out_file):
    #make dictionary with all the dataframes of different distances
    data_dict = {}
    for file in file_list:
        name = file.split("/")[-1]
        #only read in cosine and Pseudo negative set.
        data_dict[name] = pd.read_csv(file, engine = 'python', sep = "\t", index_col = False, \
            usecols=['pair', 'cosine', 'Interaction']).query('Interaction != "RusselNeg"')
        data_dict[name].loc[:,'cosine'].astype(float)
    #data_dict['egg_hmm_distances_5pubID']#.query('Interaction != "BioGrid"').sort_values(by = ['Interaction'], ascending = True)
    #confidence intervals
    c_ints = {} #confidence intervals dictionary (key = orthology, value are the confidence intervals)
    #permute_auc = {} #
    for name in data_dict:
        #default boots/perms = 1000, bootstrap the AUC values
        #since we have an extreme positive and negative set size difference,
        #do we need to make a selection of the big negative set first?
        #or is here like John said the impalance also not a problem when calculating the AUC?
        c_ints[name] = bootstrap_auc(data_dict[name].loc[:,'Interaction'],\
        data_dict[name].loc[:,'cosine'])

    for name, value in c_ints.items():
        print("%s\t%f\t%f\n" % (name, value[0], value[1]))
    return c_ints

if __name__ == '__main__':
#    Process(target=bootstrap_all, args =(file_list,)).start()
    Process(target=bootstrap_all, args =(file_list[0:2],out_file,)).start()
    Process(target=bootstrap_all, args = (file_list[2:4],out_file,)).start()
    Process(target=bootstrap_all, args = (file_list[4:6],out_file,)).start()
    Process(target=bootstrap_all, args = (file_list[6:8],out_file,)).start()
