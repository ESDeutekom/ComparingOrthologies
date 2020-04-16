#python 3
'''
    Author: John van Dam.
    Created: July 28th, 2017
    Adjusted: August 28th, 2017
    Version: 0.0.0
    This version adjust by Eva Deutekom
    Adjustments include:
    	1. Reading tab delimeted binary profile files instead of binary vector per profile
    	2. Line 152: how tree is read in: format = 1.
    	3. internal node name numbering removed, this seems to overwrite clade node names?
        Reference: Deutekom et al. 2019 PLOS Computational Biology 15(8): e1007301.
    Purpose of script: The purpose of this script is to calculate
    dollo parsimony for a set of phylogenetic profiles given a species
    tree.

    This is a partial reimplementation of Phylip's dollo parsimony script.
    Reference: Kensche et al. 2008 J R Soc Interface 5(19): 151-70.
    Note: ete node.add_feature() converts all values to string...
'''
import sys
import os
import argparse
import logging
#import re
#from time import sleep
import ete3
from tqdm import tqdm  # A progress bar, use as for i in tqdm(iterable):

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
# logger.setLevel("DEBUG")

# Functions!
def checkArgFileExists(path):
    ''' For ArgumentParser, Checks if path is actually a file. '''
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("File %s does not exist." % path)
    elif not os.access(path,os.R_OK):
        raise argparse.ArgumentTypeError("File %s is not readable. Check your permissions!" % path)
    else:
        return path

def fileExists(path):
    if os.path.exists(path):
        return True
    else:
        return False

def fileIsReadable(path):
    if fileExists(path) & os.access(path,os.R_OK):
        return True
    else:
        return False

def fileIswritable(path):
    if fileExists(path) & os.access(path,os.W_OK):
        return True
    else:
        return False

def checkOutputFile(path):
    if fileExists(path):
        # If file exists, we need to check if force is set
        if args.force:
            if fileIswritable(path):
                # This is good
                return True
            else:
                # This is not good
                raise argparse.ArgumentTypeError("File %s is not writable. Check your permissions!" % path)
        else:
            # Force is not set so throw error, clustalo style!
            logger.error("Cowardly refusing to overwrite existing file %s! Set --force to force overwriting existing files." % path)
            exit(0)
    else:
        return True

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).
    The "answer" return value is one of "yes" or "no".
    From: https://gist.github.com/hrouault/1358474
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def match_nodes(treeA,treeB):
    '''
    This function creates a list of tupules to match any node in tree A to the
    corresponding node in tree B. I am assuming identical tree topology for both
    trees here.
    '''
    return zip(list(treeA.traverse()),list(treeB.traverse()))

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("-p","--profile_file", metavar="PROFILE", required=True, type=checkArgFileExists, help="File containing the phylogenetic profiles.")
parser.add_argument("-s","--species_file", metavar="SPECIES", required=True, type=checkArgFileExists, help="File containing the order of the species in the profiles.")
parser.add_argument("-t","--tree_file", metavar="TREE", required=True, type=checkArgFileExists, help="File containing the species tree [newick].")
parser.add_argument("-po","--profile_out_file", metavar="PROFILE_OUT", type=str, default='profile_out.txt', help="Filename to write detailed profiles out to.")
parser.add_argument("--force", action='store_true', help="Force overwriting of the profile out file and output file.")
parser.add_argument("--debug", action='store_true', help="Also print debugging information.")
args = parser.parse_args()

# Check arguments
if args.debug:
    logger.setLevel("DEBUG")

## Check if the profile out file exists and warn for overwriting unless --force is set
checkOutputFile(args.profile_out_file)
## Some messages for the user concerning his/her choices.
if args.force:
    logger.warning("--force is set. The file %s will be overwritten!" % args.profile_out_file)

logger.info("Commencing...")


logger.info("Loading the list of species that are in the profiles.")
# we're assuming one species name per line.
profile_species = list()
with open(args.species_file, 'r') as sfh:
    for species in sfh:
        profile_species.append(species.rstrip())
logger.debug(profile_species)


logger.info("Loading and preparing species tree.")
base_tree = ete3.Tree(args.tree_file,format=1)  # Load the species tree.
tree_species = set(base_tree.get_leaf_names())
print(len(tree_species))
logger.debug(tree_species)
# Check if the species in the tree matches the species in the profile.
if set(profile_species) == tree_species:
    logger.info("The species in the profiles match the species in the tree.")
else:
    logger.warning("Species in the profiles do not match the species in the tree!")
    if len(set(profile_species).difference(tree_species)):
        logger.warning("Species which are in the profiles, but are not found in the tree: " + ", ".join(set(profile_species).difference(tree_species)))
        logger.warning("These species will be ommited in the analysis.")
    if len(tree_species.difference(set(profile_species))):
        logger.warning("Species which are in the tree, but are not seen in the profiles: " + ", ".join(tree_species.difference(set(profile_species))))
        base_tree.prune(tree_species.intersection(set(profile_species)))
        pruned_tree_file = args.tree_file+".prunned"
        base_tree.write(outfile=pruned_tree_file, format=9)
        logger.warning("Prunned the species tree to match the profiles and saved it to %s." % pruned_tree_file)
# Add internal node names
#node_number = 0
#for node in base_tree.traverse("levelorder"):
#    if not node.is_leaf():
#        node_number += 1
#        node.name="n%i" % node_number


logger.info("Loading profiles.")
basic_profiles = list()
with open(args.profile_file) as pfh:
    next(pfh) #skip header line
    for profile in pfh:
        # Parse the profile
        lines = profile.rstrip().split("\t") #for tab delimeted file
        profile_name = lines[0]
        profile_vector = "".join(lines[1:])
        if (len(profile_vector) != len(profile_species)):
            logger.error("The absence/presence profile for %s is not the same size as the number of species in %s!" %(profile_name,args.species_file))
            exit(1)
        profile_dict = dict(zip(profile_species,profile_vector))
        basic_profiles.append((profile_name,profile_dict))
logger.info("Loaded a total of %i profiles." % len(basic_profiles))

logger.debug("Setting up the ")

logger.info("Performing Dollo parsimony analysis.")
try:
    pofh = open(args.profile_out_file,'w')
    pofh.write("Profile_name\tExtended_Newick_tree\n")
except:
    logger.error("Unable to write to the profiles output file %s." % args.profile_out_file)
    exit(1)
for (profile_name,profile_dict) in tqdm(basic_profiles, dynamic_ncols=True):  #tqdm to create progress bar
    profile_tree = base_tree.copy(method='newick')
    # Now go through each level of the tree starting at the leaves to ancestrally recontruct the absence/presence for internal nodes
    for node in profile_tree.traverse(strategy='postorder'):  # The postorder is the order of tree traversal that is right for my purpose!
        # Check if node is a leaf node.
        logger.debug('Working on node %s' % node.name)
        if node.is_leaf():
            try:
                node.add_feature('presence',profile_dict[node.name])
            except:
                logger.error("Unable to set the presence for species '%s' in profile tree for profile '%s'" %(node.name,profile_name))
        else:
            # Node is not a leaf, so we will need to get it's children to determine the presence feature
            # Assuming both bifurcating and multifurcating tree...
            node_children_presence = list()
            for child in node.get_children():
                node_children_presence.append(child.presence)
            logger.debug("state of the children are: " + str(node_children_presence))
            #node_children_presence = [child.presence for child in node.get_children()]
            if set(['0']) == set(node_children_presence): # all children are absent, so this node is also absent
                node.add_feature('presence','0')
            elif set(['0','1']) == set(node_children_presence): # one child is present, but the other is not, so basically we can't tell yet
                if node_children_presence.count('1') >= 2:  # aka if node is multifurcating and more than one child is present.
                    node.add_feature('presence','1')
                else:
                    if node.is_root():
                        node.add_feature('presence','0')
                    else:
                        node.add_feature('presence','?')
            elif set(['1']) == set(node_children_presence):  # all children are present, so this one's easy
                node.add_feature('presence','1')
            elif set(['0','?']) == set(node_children_presence): # One child is absent, one child is undetermined, so we can not tell yet unless we are root (but we'll save this for later)
                if node_children_presence.count('?') >= 2:  # aka if node is multifurcating and more than one child is present.
                    node.add_feature('presence','1')
                else:
                    if node.is_root():
                        node.add_feature('presence','0') # Because if one child is absent and one child is present assume "gain" in the child.
                    else:
                        node.add_feature('presence','?')
            elif set(['1','?']) == set(node_children_presence): # One child is present, one child is undetermined, so we assume this node is present.
                node.add_feature('presence','1')
                # We'll have to update all children's ND's too, but we are saving that for later
            elif set(['?']) == set(node_children_presence): # Both children are ?'s meaning there are presences on both branches, so this node is also present
                node.add_feature('presence','1')
            elif set(['0','1','?']) == set(node_children_presence): # This is a special case that can only occur in multifurcating trees
                # Because we do not know the actual order of descent for this node we have some kind of conundrum.
                # To solve this we are going to assume that if "at least" two sister branches has a presence, the parental node has also a presence.
                node.add_feature('presence','1')
            else:
                logger.error("Profile %s: Unable to determine absence/presence state for node %s. This should not be possible!" % (profile_name,node.name))
                exit(1)
        logger.debug('State of node %s is set to %s' % (node.name,str(node.presence)))

    # Now we need to traverse the tree again, but in reverse order to fix the ?'s by copying the ancestor state to the current node with undefined presence/absence.
    logger.debug("Traversing the tree in oposite direction to fill in the blanks and determine events.")
    for node in profile_tree.traverse(strategy='preorder'):
        # Fill in absence and presences
        if node.presence == '?':
            if node.up.presence == "1":  # 1 as a string because node.add_feature() converts everything to strings
                node.add_feature('presence','1')
            elif node.up.presence == "0":
                node.add_feature('presence','0')
            elif node.up.presence == '?':
                logger.error("While traversing the tree I found an undetermined ancestral node. This should not be possible! Node: %s" % node.name)
                exit(1)
            else:
                logger.error("While traversing the tree I found an ancestral node without ancestral state set. This should not be possible! Node: %s, Ancestral node: %s" % (node.name, node.up.name))
                exit(1)
        # Determine events at the node
        # Because we completed all dependencies for this node to call events... lets call the events...
        # Doing it now saves another tree traversing cpu consuming step.
        if node.is_root():
            if node.presence == '1':
                node.add_feature('event','gained or inherited')  # If a gene is already present at root, assume gain (I thought this was part of the algorithm, apparently not, see below)
            else:
                node.add_feature('event','none')  # Gene is not invented yet
        else:
            if node.presence == '1':
                if node.up.presence == '1':
                    node.add_feature('event','inherited')  # vertical inheritance
                else:
                    node.add_feature('event','gained') # First it wasn't there, but now it is! Gain!
            else:
                if node.up.presence == '1':
                    node.add_feature('event','loss') # First is was there, but now it isn't. Loss!
                else:
                    if (node.up.event == 'loss') | (node.up.event == 'ancestral_loss'):
                        node.add_feature('event','ancestral_loss')
                    elif node.up.event == 'none':
                        node.add_feature('event','none')  #gene has not been invented yet
                    else:
                        logger.error("Combination of Ancestral events and current state not accepted. Something went wrong!")
                        exit(0)
        if node.is_leaf():
            logger.debug("%s: %s state is %s by event %s" % (profile_name,node.name,node.presence,node.event))
        if node.is_root() & (node.presence == '1'):
            logger.debug("%s is a LECA OG!" % profile_name)
    # Printing profile name and presence absence acestral reconstruction in newick format
    logger.debug(profile_name + "\t" + profile_tree.write(features=["presence",'event'],format=3,format_root_node=True))
    logger.debug(profile_tree.get_ascii(attributes=('presence','event')))
    pofh.write("%s\t%s\n" % (profile_name,profile_tree.write(features=["presence",'event'],format=3,format_root_node=True)))
pofh.close()




logger.info("All done!")
# The end!
