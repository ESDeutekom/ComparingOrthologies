#python3

#####################################################################################
## Remove empty columns from (MSA) files left over by MSA containing more sequences
## that were parsed out for stricter LECA definition
#####################################################################################
import sys
from Bio import AlignIO
import re
import sys

if len(sys.argv) != 3:
    print("Need 2 arguments: [MSA file] [MSA trimmed out file name]")
    sys.exit()

msa = sys.argv[1] #old MSA files
msa_out = sys.argv[2] #new MSA files

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()


#alignment object
alignment = AlignIO.read(msa, "fasta")
#Columns to keep
col_keep_num = [] #columns list that are not empty
for col in range(alignment.get_alignment_length()):
    if re.search('[a-zA-Z]', str(alignment[:,col])): #anything that has aminoacids
        col_keep_num += [col] #keep columns

msa_out = open(msa_out, "w")
for s in alignment:
    msa_out.write(">%s\n" % s.id) #sequence ID
    msa_out.write("%s\n" % "".join([s.seq[i] for i in col_keep_num])) #Sequence with columns to keep
msa_out.close()
