#python3

################################################################################
## Getting PANTHER LECA's from ancestral genomes trees by Huang et al. 2018
## Using additional parsing criteria to be more strict
################################################################################
import sys
from ete3 import Tree

if len(sys.argv) != 5:
    print("Need 4 arguments: [ancestral genome tree] [out file directory name] [Patnher ID of tree] [stats out file name]")
    sys.exit()

tree_file = sys.argv[1] #ancestral genome trees
out_file = sys.argv[2] #outfile directory
pthrID = sys.argv[3] #Ancestral panther ID
stats_file = sys.argv[4] #some stats on the ogs

try:
	open(sys.argv[1])
except IOError:
    print("No such input file"); sys.exit()


supergroupsL = ["Excavates","Viridiplantae", "Alveolata-Stramenopiles" ]
supergroupsR = ["Amoebozoa", "Opisthokonts"]

#tree file has additional lines at the end containing info of leafs
#only need first line that contain tree line
#if leafs have -1 branch length, this means they are lost and do not contain data
#(i.e. empty sequences in panther msa files)
tree = open(tree_file, "r").readline()

inRoot = "unknown"
t = Tree(tree)
#criteria for in LECA, left and right, in atleast 3 supergroups in total and in the root
leaf_listL = []
leaf_listR = []
leca_dict = {}
for node in t.traverse():
    if hasattr(node,"S"): #NHX file that has attributes
        name = node.S
        ID = node.ID
        if name ==  "Eukaryota":
            LecaID = ID
            leaf_listL = []
            leaf_listR = []
            child_listL = []
            child_listR = []
            for node_child in node.traverse():
                if hasattr(node_child,"S"):
                    name_child = node_child.S
                    ID_child = node_child.ID
                    if name_child in supergroupsL:
                        for leaf in node_child: #all descendants from this node are LECA
                            if float(leaf.dist) != -1:
                                leaf_listL += [leaf.name]
                        if leaf_listL:
                            if name_child not in child_listL:
                                child_listL += [name_child]
                    elif name_child in supergroupsR:
                        for leaf in node_child: #all descendants from this node are LECA
                            if float(leaf.dist) != -1:
                                leaf_listR += [leaf.name]
                        if leaf_listR:
                            if name_child not in child_listR:
                                child_listR += [name_child]
            if all([child_listL, child_listR]): # if left and right of the root/eukaryota node
                if (len(child_listL) + len(child_listR) >= 3): # if in atleast 3 supergroups
                    if (len(leaf_listL) + len(leaf_listR) >= 3): # if has atleast 3 sequences
                        file_name = out_file + pthrID + "_" + LecaID #file name containing Panther fam ID and LECA ID
                        file_out = open(file_name, "w")
                        file_out.write("\n".join(["\n".join(leaf_listL), "\n".join(leaf_listR)]) + "\n")
                        #Some stats on the panther leca's
                        stats = open(stats_file, "a")
                        stats.write("\t".join([pthrID, LecaID, ",".join(child_listL), str(len(leaf_listL)), ",".join(child_listR), str(len(leaf_listR))]) + "\n")
