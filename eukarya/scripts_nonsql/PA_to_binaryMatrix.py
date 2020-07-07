#python3

import os
import sys
import pandas as pd

################################################################
## Make a binary matrix from the Presence and Absence matrix
###############################################################

if len(sys.argv) != 3:
    print("Need 2 arguments: [Presence and absence input file] [binary output file name]")
    sys.exit()

PA_file = sys.argv[1]
binary_outfile = sys.argv[2]

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()

#Check if file is not empty
if os.path.getsize(sys.argv[1]) <= 1:
    print(PA_file, "file is empty"); sys.exit()

#read in the PA files
df = pd.read_table(PA_file, index_col = 0, sep = "\t", header = None, names = range(0,170), engine = "python") #read in dataframe from PA file
#quite memory intensive!
binary_matrix = pd.get_dummies(df.unstack()).sum(level=1) #unstack to pivot data, sum

binary_matrix.to_csv(binary_outfile, sep = "\t")
