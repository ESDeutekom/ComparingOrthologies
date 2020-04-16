#python3

import sys
from ete3 import Tree

########################################################################################
## Find LECA orthologous groups
## This is the strict definition of LECA so presense must be in atleast 3 supergroups
## Atleast in the root, left and right of the root (aka unikont / bikont)
## Hardcoded to trees in eukarya
########################################################################################

if len(sys.argv) != 2:
    print("Need 1 arguments: [dollo tree file from dollo parsimony]")
    sys.exit()

dollo_singlefile = sys.argv[1]

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()

#define the supergroups
#Because GTHE and TTRA do not have (cannot have) a node annotation
#they are added seperatly
supergroupsLeft = ["Excavata","Archaeplastida","SAR","Haptophyta","GTHE"]
supergroupsRight = ["Amoebozoa","Opisthokonta","TTRA"]
#read in tree
dollo_file = open(dollo_singlefile, "r")
dollo_file.readline()
for line in dollo_file:
    line = line.rstrip().split('\t')
    OG_id = line[0]
    treeLine = line[1]
    tree = Tree(treeLine, format = 1)
    # for tree, get the node and check if root has a presence 1
    #additionally, see if branches/nodes are in supergroups and see if presence is at least in 3 supergroups
    superCountleft = 0 # counter for the left side of the root to check if there is a left side
    superCountright = 0 # counter for right side
    rootFind = "unknown" # See if there is a root

    for node in tree.traverse(strategy='preorder'):
        if (node.name in supergroupsLeft) & (node.presence == '1'):
            superCountleft += 1 # if there is a presence in left root --> +1
        elif node.is_root() & (node.presence == '1'):
            rootFind = "yes" # if there is a root, rootFind --> yes
        elif (node.name in supergroupsRight) & (node.presence == '1'):
            superCountright += 1 # if there is a presence in right root --> +1
    superCount = superCountright + superCountleft # total supergroups must be 3 or higher
    #if root is found, presence is in both left and right root, and total supergroup presence is 3 or largers, we found a LECA PFAM
    if (rootFind == "yes") & (superCount >= 3) & (superCountleft >= 1) & (superCountright >= 1):
        print(OG_id)
