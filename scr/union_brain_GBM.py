#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:06:45 2022

Compare GEM (brain and tumor) and do a UNION

@author: johanna
"""
import os
import pandas as pd


def ordertup_genemet(lot):
    """
    input: list of lists : [['asparagine|C_s', 'SLC7A6'], ['ABCA',... ],...]

    output : list of tuples : [('SLC7A6', 'asparagine|C_s') , ... ]
    """
    outtl = list()
    for l in lot:
        if "|" in l[0]:
            outtl.append((l[1],l[0]))
        else:
            outtl.append(tuple(l))
    return outtl


def messagedifftwoEdgesSets(storeEds, index0, index1):
    # index0, index1 : the two elements of the dictionary to compare
    name1 = list(storeEds.keys())[index0]
    name2 = list(storeEds.keys())[index1]
    ostr = ''
    diff12 = set(storeEds[name1])  - set(storeEds[name2])    
    ostr += f"Edges present in {name1}, but not in {name2}\n"
    ostr += f"{name1} - {name2} = {len(diff12)}\n"
    for i in list(diff12):
        ostr += "\t".join(i) +"\n"
    ostr += f"[END {name1} - {name2}]\n\n"
    diff21 = set(storeEds[name2]) - set(storeEds[name1])
    ostr += f"Edges present in {name2}, but not in {name1}\n"
    ostr += f"{name2} - {name1} = {len(diff21)}\n"
    for i in list(diff21):
        ostr += "\t".join(i) +"\n"
    ostr += f"[END {name2} - {name1}]"
    
    return ostr
    
### START
print(" ** GEMs : edges comparison and Union ** ")
# ** took trimmed and did a UNION first
dirgems = os.path.expanduser("~/recast_GEMpub/")
ffix = "union_brain_GBM"
os.chdir(dirgems)
existinggems = ['brain_normal', 'GBM_tumor']
ofirep = f"{ffix}/report_{existinggems[0]}_{existinggems[1]}.txt"
ofileE = f"{ffix}/union_{existinggems[0]}_{existinggems[1]}.txt"
if not os.path.exists(os.path.join(os.getcwd(),ffix)):
    os.makedirs(ffix)
    
 
storeEds = dict()   

print("opening available GEMs, ordering edges (gene,met), removing duplicates")    
for i in existinggems:
    tmpf = [j for j in os.listdir(dirgems+i) if i+"_trimmed.txt" == j]
    #tmpf = [j for j in os.listdir(dirgems+i) if "bipartite_"+i+".txt" == j]
    assert len(tmpf) == 1, "ambiguous _trimmed.txt files"
    print(tmpf[0])
    with open(dirgems + i + "/" + tmpf[0], "r") as f:
        po = f.readlines()
    lol = [i.strip().split("\t") for i in po]
    #print(lol)
    # reorder pairs, make sure this order: (gene, metabolite)
    loto_ = ordertup_genemet(lol)
    # get rid of duplicates if any : set
    tmp = set(loto_); del(loto_)
    loto_ = list(tmp)
    # store edges list into a dictionnary
    storeEds[i] = loto_
     
mess1221 = messagedifftwoEdgesSets(storeEds, 0, 1)

# union 
edunion = set(storeEds['brain_normal']).union(set(storeEds['GBM_tumor']))

print("Edges GBM : ", len(storeEds['GBM_tumor']))
print("edges normal U GBM : " , len(edunion))

mess1221 += "\n============================================\n"
mess1221 += f"Edges brain (normal) : {len(storeEds['brain_normal'])}\n"
mess1221 += f"Edges GBM : {len(storeEds['GBM_tumor'])}\n"
mess1221 += f"Edges normal U GBM : {len(edunion)}\n"
with open(ofirep, "w") as f:
    f.write(mess1221)
    
uniontx = ''
for k in ["\t".join([i,j]) for (i,j) in edunion]:
    uniontx += k + '\n'

with open(ofileE, "w") as f:
    f.write(uniontx)    
    