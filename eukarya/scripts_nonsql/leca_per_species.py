#python3

import os
import sys
import random
import pandas as pd
####################################################################
## Count LECA OGs and fractions per species
####################################################################

if len(sys.argv) != 4:
    print("Need 3 arguments: [profile] [LECA list file] [out file name]")
    sys.exit()

pro_file = sys.argv[1] #bash array with file names
LECA_file = sys.argv[2] #bash array with file names
out_file = sys.argv[3] #out file

#Check if file is not empty
for file in (sys.argv[1], sys.argv[2], sys.argv[3]):
    if os.path.getsize(file) <= 1:
        print(file, "file is empty"); sys.exit()

orthology = str(pro_file).split("/")[-1].split("_binary")[0]

#read in leca OGs
LECA_file = open(LECA_file, "r")
leca_list = []
for line in LECA_file:
    leca = str(line.rstrip())
    leca_list += [str(leca)]
leca_size = len(leca_list)

#start counting the presences and absences per species for the leca OGs
pro_file_df = pd.read_csv(pro_file, sep = "\t", index_col = 0)
og_size = len(pro_file_df.index)
pro_file_df.index = pro_file_df.index.astype(str)

pro_file_df_leca = pro_file_df[pro_file_df.index.isin(leca_list)].sum(axis = 0) #only select leca OGs and sum
summed_df_leca = pro_file_df_leca.rename(orthology+"_leca_ogs_of_"+str(leca_size)) #give it the orthology name


pro_file_df_nonleca = pro_file_df[~pro_file_df.index.isin(leca_list)] #only select ones not in leca
#need to remove not group OGs (singlets are still in here)
print(len(pro_file_df_nonleca))
summed_df_ogs = pro_file_df_nonleca.sum(axis = 1) #Sum over the OGs to see if they contain more that 2 species

summed_df_ogs_names = list(summed_df_ogs[summed_df_ogs > 1].index) #get the OG names
print(len(summed_df_ogs_names))

#summed_df_ogs_namesw = list(summed_df_ogs[summed_df_ogs == 1].index)
#summed_df_ogs.columns = summed_df_ogs.columns.astype(str)
summed_df = pro_file_df_nonleca[pro_file_df_nonleca.index.isin(summed_df_ogs_names)].sum(axis = 0)#sum over all species the precenses in the OGs
print(summed_df)
#summed_df.columns = summed_df.columns.astype(str)
summed_df = summed_df.rename(orthology+"_ogs_of_"+str(og_size))

df = pd.read_csv(out_file, sep = ",", index_col = 0)
df = pd.concat([df, summed_df_leca], axis=1)#add this to already existing files/calculations
df = pd.concat([df, summed_df], axis=1)
df.to_csv(out_file)#, header = False)
