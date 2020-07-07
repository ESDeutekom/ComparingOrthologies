#python3

import os
import sys
import pandas as pd
import numpy as np
import re
from itertools import combinations
import copy as cp
best_OGs = sys.argv[1]

df = pd.read_csv(best_OGs, sep = ",", index_col = 0)
df_dict = {}
orth_list = []
for column, rows in df.iteritems():
    column = re.split(", |'|\(|\)", column)
    column_name = [el for el in column if el]
    #print(column_name)
    orth1 = column_name[0]
    orth2 = column_name[1]
    for row in rows:
        if isinstance(row, str):
            row = re.split(", |'|\(|\)", row)
            row_name = [el for el in row if el]
            og1 = row_name[0] +"."+ orth1
            og2 = row_name[1] +"."+ orth2
            #print(orth1, orth2, og1, og2)
            if orth1 not in df_dict:
                df_dict[orth1] = {}
                df_dict[orth1][og1] = [og1] #add self
                df_dict[orth1][og1] += [og2]
            else:
                if og1 not in df_dict[orth1]:
                    df_dict[orth1][og1] = [og1] #add self
                    df_dict[orth1][og1] += [og2]
                else:
                    if og2 not in df_dict[orth1][og1]:
                        df_dict[orth1][og1] += [og2]
            if orth2 not in df_dict:
                df_dict[orth2] = {}
                df_dict[orth2][og2] = [og2] #add self
                df_dict[orth2][og2] += [og1]
            else:
                if og2 not in df_dict[orth2]:
                    df_dict[orth2][og2] = [og2] #add self
                    df_dict[orth2][og2] += [og1]
                else:
                    if og1 not in df_dict[orth2][og2]:
                        df_dict[orth2][og2] += [og1]

def make_single_dict(nested_dict, key):
    new_dict = nested_dict[key]
    print(key)
    return new_dict

ob = make_single_dict(df_dict, 'orthofinder_blast_e-3')
od = make_single_dict(df_dict, 'Swiftortho_c50')

df_copy = cp.deepcopy(df_dict)

new_compared_dict = {}
orthology_list = []
for orthology in df_dict:
    orthology_list += [orthology]
    #combi_list = list(combinations(orthology_list, 2))


#everything the same must be removed out of ortholgy dict and added once to over new compared dict
#The ones left in the ort
groups_list = []
for orthology in orthology_list:
    x = list(df_dict[orthology].values())
    for groups in x:
        groups = set(groups)
        if groups not in groups_list:
            groups_list += [groups]
            print(groups_list)

    #new_compared_dict = {k: x[k] for k in x if k in y and x[k] == y[k]}
#print(new_compared_dict)




"""


for i in range(1,len(orth_list))[::-1]: #the first has all comparisons, last is not it's own dict
    ortho_n = orth_list[i]
    ortho_m1 = orth_list[i-1]
    print(ortho_n, ortho_m1)
    #print(df_dict[ortho_n])

    #if ortho in df_dict:

    #    print(df_dict[orth_list[0]].values() > 2)
        #print(df_dict[ortho].values())
        #intersect_org = df_dict[orth_list[0].values()].intersect(df_dict[ortho].values())
        #if  len(intersect_org) > 1:
        #    print(intersect_org)

ortho_all = {}

for ortho1 in df_dict:
    ortho_all[ortho1] = {}
    for ortho2 in df_dict[ortho1]:
        ortho_all[ortho1][ortho2] = set(df_dict[ortho1][ortho2].keys())
        #ogs2 = [list(el.keys())[0] for el in df_dict[ortho1][ortho2].values()]

for ortho, tuple_dict in ortho_all.items():
    for i in range(2,len(tuple_dict.keys())+1):
        combis = list(combinations(tuple_dict.keys(), i))
        for combi in combis:
            list_list = []
            print(ortho, combi)
            for el in combi:
                list_list += [tuple_dict[el]]
            print(len(list(set.intersection(*map(set, list_list)))))
            #print(len(ortho_all[ortho][el[0]] & ortho_all[ortho][el[1]]))
            #print(len(ortho_all[ortho][el[0]].intersection(ortho_all[ortho][el[1]])))
    #for og_list in value:
"""
