#python3
#################################################################################
## Parse out empty sequences from MSA files Ancestral genomes Huang et al. 2018
## Parse out non LECA/Ancestral Panther IDs
#################################################################################
import sys
from Bio import SeqIO
import re

if len(sys.argv) != 3:
    print("Need 2 arguments: [MSA file] [Panther LECA IDs]")
    sys.exit()

msa = sys.argv[1] #msa file to parse for leca panthers and remove empty/lost sequences
leca_id = sys.argv[2] #LECA panther ID

try:
	open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

with open(msa, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id == leca_id: #if record.id is LECA panther id
            seq = record.seq
            if re.search('[a-zA-Z]', str(seq)) == None: #no characters means empty sequence string
                pass
            else:
                print(">"+record.id+"\n"+seq ) #print out to file with panther + leca id in wrapper
#Example wrapper
#for dir in  ./PantherLECAslists/*
#do
#    lecaID=$(basename -- "$dir"); pthrID=${lecaID%_*}
#    echo $lecaID, $pthrID
#    while read LINE
#    do
#        echo $LINE
#        ./parse_nonSeq.py ../$pthrID/tree.mia $LINE >> ./PantherLECA_MSA/$lecaID
#    done < ./PantherLECAslists/$lecaID
#done
