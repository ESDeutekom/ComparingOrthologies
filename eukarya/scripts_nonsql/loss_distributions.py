#python 3

import os
import sys
from find_loss import *

if len(sys.argv) != 4:
    print("Need 3 arguments: [Dollo parsimony tree] [LECA OG list] [outfile name]")
    sys.exit()

dollo_tree = sys.argv[1]
leca_file = sys.argv[2]
out_file = sys.argv[3] #out file"""

try:
    open(sys.argv[1])
    open(sys.argv[2])
except IOError:
    print("No such input file"); sys.exit()

_,indep_loss = loss_dict(dollo_tree, leca_file)

loss_out = open(out_file, "w")
for key, value in indep_loss.items():
    loss_out.write('%s\t%.0f\n' % (key, value))
loss_out.close()
