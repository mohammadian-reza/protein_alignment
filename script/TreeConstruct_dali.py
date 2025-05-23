# conda install -c bioconda fastme
# pip install ete3 colour wget statsmodels matplotlib seaborn

import glob
import ete3
import wget 
import pandas as pd
import os, sys
import time
import json
from src import AFDB_tools, treescore , foldseek2tree
import colour
import numpy as np
import toytree
import pickle
from matplotlib import pyplot as plt
from Bio import SeqIO
import toytree
import toyplot.svg
import csv

method = 'dali' # deepblast, TMalign, deepalign
sim = 'zscore' # or 'TM-score'

cluster = ''

draw = False
#get trees
swisstreepath = f'../SwissTree{cluster}/*/*nhxpruned.nhx'
swiss_trees = glob.glob(swisstreepath)

import subprocess, shlex
def madroot(t):
    args = 'madroot/mad ' + t
    subprocess.run(shlex.split(args))
    return t+ '.rooted'

## calculate distance
import traceback
import seaborn

treescores_st = []
treescores_str = []
treescores_struct = []
rfdistances = []

overwrite = False

if not os.path.exists(f'../SwissTree{cluster}.congruence') or overwrite:
    with open(f'../SwissTree{cluster}.congruence', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'Congruence_score'])

if not os.path.exists(f'../SwissTree{cluster}.rf') or overwrite:
    with open(f'../SwissTree{cluster}.rf', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'RF_distance','max_rf'])

#for t in swiss_trees:
for t in swiss_trees:
    folder = os.path.dirname(t)
    swisst = pickle.load(open(folder+'/swiss_toytree.pkl','rb'))
    out_tree_kernel = folder + '/dali_unrooted.nhx'
    if not os.path.exists(out_tree_kernel):
        print(f'Tree not found in {folder}')
        continue
    os.system(f'ls {folder}/structs | sort | uniq > {folder}/identifiers_order.txt')
    tre=toytree.tree(out_tree_kernel,tree_format=0)
    with open(f'{folder}/identifiers_order.txt','r') as f:
        ids = [i.strip().replace('.pdb','') for i in f]

    name_dict={}
    for i,n in enumerate(ids):
        name_dict['s'+str(i+1).zfill(3)+'A']=n
    tre=tre.set_node_values(feature='name',values=name_dict)
    tre = ete3.Tree(tre.write())
    with open(folder+'/dali_unrooted_rename.nhx','w') as treeout:
        treeout.write(tre.write(format=1))
    out_tree_kernel = folder+'/dali_unrooted_rename.nhx'

    out_tree_kernel = foldseek2tree.postprocess(out_tree_kernel, out_tree_kernel.replace('.nhx','_PP.nhx'))
    out_tree_kernel = madroot(out_tree_kernel)
    tre = toytree.tree(out_tree_kernel)
    swiss_tre = madroot(t)
    swiss_tre = toytree.tree(t)
    
    unidf = AFDB_tools.grab_entries(ids)
    lineages = treescore.make_lineages(unidf)
    
    tre = treescore.label_leaves(tre,lineages)
    swiss_tre = treescore.label_leaves(swiss_tre,lineages)
    overlap = treescore.getTaxOverlap(tre.treenode)
    
    try:
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode )
        rfdistances.append(rf)
    
    except:
        print('rferr')
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode,unrooted_trees=True )
        rfdistances.append(rf)
            
    print('struct score: ', tre.treenode.score, 'ST score: ' , swisst['tree'].treenode.score )

    treescores_struct.append( tre.treenode.score)
    treescores_st.append( swisst['tree'].treenode.score  )

    
    with open(f'../SwissTree{cluster}.congruence', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], f'{method}_{sim}', len(tre), tre.treenode.score])

    with open(f'../SwissTree{cluster}.rf', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], f'{method}_{sim}', len(tre), rf[0], rf[1] ])
    
if draw:
    plt.scatter( x = treescores_struct , y = treescores_st  , c = 'b', ls = '-')
    plt.plot( [0,2000] , [0,2000])
    plt.xlabel('structure')
    plt.ylabel('swisstree')
    plt.show()

