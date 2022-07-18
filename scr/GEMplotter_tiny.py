#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:49:20 2022

@author: johanna
"""
import pickle
import os
import argparse
import networkx as nx
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")


#os.chdir(args.mywdir)

mywdir = "~/recast_GEMpub/"
os.chdir(os.path.expanduser(mywdir))

with open("GBM_tumor/dictioFull_GBM_tumor.pkl" , "rb") as f:
    tumorD = pickle.load(f)
    
erp = tumorD['enzy_reac_prod']
erp[6]

G = nx.Graph()


for k in range(0,1000):
    edgeshere = []
    for z in erp[k]['enzymes']:
        for y in erp[k]['reactants']:
            for x in erp[k]['products']:
                G.add_node(z, atype='enzyme', shape='s', color='green')
                G.add_node(y, atype='metabolite',shape='o' , color='orange')
                G.add_node(x, atype='metabolite', shape='o', color='orange')
                edgeshere.append((z,y)) # enzyme with reactant
                edgeshere.append((z,x)) # enzyme with product
    print(len(edgeshere))
    edgescleared = list(set(edgeshere))  # remove repeated edges
    print(edgescleared)
    G.add_edges_from(edgescleared)

enzymes = [n for (n, aty) in nx.get_node_attributes(G,'atype').items() if aty == "enzyme"]
mets = [n for (n, aty) in nx.get_node_attributes(G,'atype').items() if aty == "metabolite"]
fig, ax = plt.subplots(figsize = (35,35))
#plt.figure(figsize = (35,35))
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, nodelist=enzymes, node_color='green', 
                       alpha=0.5,node_shape = 's')
nx.draw_networkx_nodes(G, pos, nodelist=mets, node_color="orange",
                       alpha=0.5, node_shape = 'o')
nx.draw_networkx_edges(G, pos, edge_color="white", alpha=0.7)
ax.set_facecolor('lightgray')
plt.savefig('yeeee', dpi=300)

enzymes
