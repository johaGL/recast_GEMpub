#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:25:08 2022

@author: johanna
"""

import os
import pickle
import argparse
import networkx as nx
from networkx.algorithms import bipartite

def lists2edges(list1, list2):
    oli = []
    for i in list1:
        for j in list2:
            if [i,j] not in oli:
                    oli.append([i,j])
    return oli

def allmetaboedges(gem):
    edges_ = []
    for k in gem['enzy_reac_prod'].keys():
        #print(f'parsing reaction {k}')
        enz_ = gem['enzy_reac_prod'][k]["enzymes"]
        rac_ = gem['enzy_reac_prod'][k]["reactants"]
        prod_ = gem['enzy_reac_prod'][k]["products"]
        subedges1 = lists2edges(enz_, rac_)
        subedges2 = lists2edges(enz_, prod_)
        edges_ += subedges1
        edges_ += subedges2
    return edges_

def doreportbipartite(suffix, GEM, edges_):     
     A = f"---------------------------------------- \
            \n{suffix}\n----------------------------------------"
     a = "\nTransformed GEM dictionnary into list of Edges\n"
     b = f"number of reactions : {len(GEM['enzy_reac_prod'])}\n"
     c = f"Edges : {len(edges_)}\
       \nNodes\
       \n\tgenes : {len(GEM['genes'])}\n\tmetabolites : {len(GEM['metabolites'])}"
     return [A, a, b, c] 
 
def giveatoy():
    pretoyD = { 4660: {'enzymes': ['SLC35D1'], 
                     'reactants': ['UDP-glucuronate|C_c', 'UMP|C_r'], 
                     'products': ['UDP-glucuronate|C_r', 'UMP|C_c']}, 
             4661: {'enzymes': ['SLC35D2'], 
                    'reactants': ['UDP-N-acetylglucosamine|C_c', 'UMP|C_r'], 
                    'products': ['UDP-N-acetylglucosamine|C_r', 'UMP|C_c']}, 
             4662: {'enzymes': ['SLC19A1'], 'reactants': ['formate|C_c'], 
                    'products': ['formate|C_r']}, 
             4663: {'enzymes': ['SLC2A1'], 
                    'reactants': ['glucose|C_c'], 'products': ['glucose|C_r']}, 
             4664: {'enzymes': ['SLC35A2'], 
                    'reactants': ['UDP-galactose|C_c', 'UMP|C_r'], 'products': ['UDP-galactose|C_r', 'UMP|C_c']}, 4665: {'enzymes': ['SLC37A4'], 'reactants': ['glucose-6-phosphate|C_c', 'Pi|C_r'], 'products': ['glucose-6-phosphate|C_r', 'Pi|C_c']}, 4666: {'enzymes': ['SLCO4A1'], 'reactants': ['estradiol-17beta 3-glucuronide|C_c'], 'products': ['estradiol-17beta 3-glucuronide|C_r']}, 4667: {'enzymes': ['AQP6', 'SLC5A1', 'AQP8', 'AQP10', 'AQP5', 'AQP7', 'AQP3', 'AQP2', 'AQP4', 'AQP1'], 'reactants': ['H2O|C_c'], 'products': ['H2O|C_s']}}
    toyD = {'enzy_reac_prod' : pretoyD}
    return toyD

    

if __name__ == "__main__":    
    
    mywdir = os.path.expanduser("~/recast_GEMpub/")
    os.chdir(mywdir)
    suffix = "GBM_tumor"
    gemfi=f"{suffix}/HMAdictio_{suffix}.pkl"
    with open(gemfi, "rb") as f:
        GEM = pickle.load(f)  
    print(GEM.keys())
    #edges_ = allmetaboedges(GEM)
    edges_ = allmetaboedges(giveatoy())

    #print(len(edges_))  
    #print(len(GEM["genes"]))
    #print(len(GEM["metabolites"]))
   
    
    sms = doreportbipartite(suffix,GEM, edges_)
    with open(f"reportbipartite_{suffix}.txt", "w") as f:
        f.writelines(sms)
    print("saved")
    G = nx.Graph()
    G.add_edges_from(edges_)
    boolbipartite = nx.is_bipartite(G)
    print(f"The graph from these Edges, is it bipartite?: *{boolbipartite}*")

    X, Y = bipartite.sets(G)
    pos = dict()
    pos.update( (n, (1,i)) for i, n in enumerate(X) )
    pos.update( (n, (2,i)) for i, n in enumerate(Y) )
    nx.draw(G, pos=pos)
    #plt.show()
    print(GEM['enzy_reac_prod'][300])
    