#python 3

import sys
from ete3 import Tree
import statistics as st
import pandas as pd
#################################################################
## Counting loss in dollo set by counting all losses in nodes
## This code does not take care of later inherited/gained genes, aka non LECA
################################################################

def loss_dict (dollo_tree, leca_file):
    leca_dict = {}
    leca_file = open(leca_file, "r")
    for lines in leca_file:
        line = lines.rstrip()
        leca_dict[line] = 0

    tree_dict = {} #dolle tree, one line per tree per orthology. irst line
    dollo_tree = open(dollo_tree, "r")
    dollo_tree.readline()
    for lines in dollo_tree:
        line = lines.rstrip().split("\t")
        OG = line[0]
        tree = line[1]
        if OG in leca_dict:
            tree_dict[OG] = tree

    for OG, tree in tree_dict.items():
        tree_structure = Tree(tree, format = 1)
        #lineage and clade loss --> independent loss
        lossClade = 0
        lossLineage = 0
        for node in tree_structure.traverse(strategy='preorder'):
            if not node.is_leaf(): #nodes that are leaf losses are not a clade. In this case node is the leaf
                if (node.event == "loss"): #if loss in node --> clade loss
                    lossClade += 1
        for leaf in tree_structure.iter_leaves():
            if (leaf.event == 'loss'): #if loss in leaf --> lineage specific loss
                lossLineage += 1
        indep_loss = lossClade + lossLineage
        leca_dict[OG] += indep_loss #add to leca dict the amount of loss for leca og

    total_indep_loss = 0
    list_val = []
    for key, value in leca_dict.items():
        list_val += [value]
        total_indep_loss += float(value)

    loss_dict = {}
    loss_dict["Independent Loss LECA OGs"]= total_indep_loss
    loss_dict["Mean independent loss"]= round(st.mean(list_val), 1)
    loss_dict["Median independent loss"] = st.median(list_val)
    loss_dict["Stdev independent loss"] = round(st.stdev(list_val),1)
    return loss_dict, leca_dict
