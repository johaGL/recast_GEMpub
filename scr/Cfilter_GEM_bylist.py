#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 12:23:37 2022

@author: johanna

examples :     
    python3 filter_GEM_bylist.py \
          --mywdir ~/recast_GEMpub/ \
          --ifile brain_normal/bipartite_brain_normal.txt \
          --odir brain_normal/ \
          --suffix brain_normal
          
   python3 filter_GEM_bylist.py \
       --mywdir ~/recast_GEMpub/ \
       --ifile GBM_tumor/bipartite_GBM_tumor.txt \
       --odir GBM_tumor/ \
       --suffix GBM_tumor
"""

import sys
import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns



parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")
parser.add_argument("--ifile")
parser.add_argument("--suffix")
parser.add_argument("--odir")
args = parser.parse_args()


#mywdir = "~/recast_GEMpub/"
os.chdir(os.path.expanduser(args.mywdir))
print(args.suffix)

#args.ifile = "brain_normal/bipartite_brain_normal.txt"
#suffix = "brain_normal"

dfi = pd.read_csv(args.ifile, header=None, sep="\t")
#print(dfi.shape)

G = nx.Graph()

for i, ro in dfi.iterrows():
    G.add_edge(ro[0], ro[1])
 
print(f'initial number of vertices: {len(list(G.nodes))}')    
degrees_ = list(G.degree())
df = pd.DataFrame(degrees_, columns = ['vertex', 'degree'])

# # degree nodes having glutamate
# for i, ro in df.iterrows():
#     if "glutamate|" in ro['vertex']:
#         print(ro['degree'])

fig, ax1 = plt.subplots()   
# see : https://stackoverflow.com/questions/33323432/add-kde-on-to-a-histogram
sns.kdeplot(ax = ax1, 
            data = df,
            color = "gray", fill=True,
            x = "degree" , 
           # log_scale = [True, False],
            alpha = 0.2)  
#ax1.set_xlim(tmpdf["ratio_value"].min(), tmpdf["ratio_value"].max())
ax2 = ax1.twinx()
sns.histplot(data= df, ax= ax2,
              x = "degree", # hue = 
              #log_scale=[True,False],
              linewidth = 0, 
              bins=200, alpha = 0.8)  
ax1.set_xlabel("vertex degree")
ax1.set_title(f'{args.suffix}')
plt.savefig(f"{args.odir}/degrees_{args.suffix}.png")



# another method:
unaccep_degree = 0
print(f"\n ** TRIMMING steps ** : \n\
      a) Using a list of promiscuous vertices, and then \n\
      b) remove  vertices with degree <= {unaccep_degree}\n")

promiscuous = ["H2O", "H+", "H2O2", "O2", "Na+", "CO2" , "CO",
               "[protein]", "phosphoprotein", "Pi", "PPi",
               "AMP", "ADP" , "ATP",
               "CMP", "CDP",  "GMP", "GDP", "UMP", "UDP",
               "NAD+","NADH", "NADP+", "NADPH", 
                "CoA", "acetyl-CoA" ] 

print(f"List of promiscuous vertices {promiscuous} \n ")

print(f"\n Step a) ")   
Gm = G.copy()
preexistingnodes = list(Gm.nodes)
promiscuous_wicompartment = set([ n for n in preexistingnodes 
                                 if n.split("|")[0] in promiscuous])
print(" nodes deleted by list", list(promiscuous_wicompartment)[:5], "...",
     list(promiscuous_wicompartment)[-5:] )

Gm.remove_nodes_from(promiscuous_wicompartment)
# now remove vertices having degree <= unaccep_degree
suppressedall = set()
condicount = 1 # random number not zero
while condicount > 0:
    degrees_m_ = list(Gm.degree())
    df_m = pd.DataFrame(degrees_m_, columns= ['vertex', 'degree'])
    df_m_nogood = df_m.loc[df_m['degree'] <= unaccep_degree]
    vx_nogood = list(df_m_nogood['vertex'])
    Gm.remove_nodes_from(vx_nogood)
    suppressed = set([i.split("|")[0] for i in vx_nogood])
    suppressedall.update(suppressed)
    condicount = df_m_nogood.shape[0]
    #print(condicount)
    #print(df_m_nogood)

print(f'\n Step b) : removing nodes with degree <= {unaccep_degree} :\n\
      {suppressedall} \n TOTAL removed : {len(suppressedall)}\n')
        
degrees_m_ = list(Gm.degree())
type(degrees_m_)
df_m = pd.DataFrame(degrees_m_, columns = ['vertex', 'degree'])


print("Verifying remaining degrees")
print(df_m.sort_values("degree", axis=0, ascending=False))
print(f"final number of vertices: {df_m.shape[0]}")


# save
df_final = list(Gm.edges())
df_final = pd.DataFrame(list(Gm.edges()), columns = ['v1', 'v2'])
df_final.to_csv(f"{args.odir}/{args.suffix}_trimmed.txt", 
                index = False, header = False, sep= "\t")

# # viewing nodes neighbouring slc16 family members
# for k in Gm.nodes:
#     if k.startswith("SLC16"):
#         print(f" -- {k}")
#         print([i for i in Gm.neighbors(k)])

# otxt = ''
# for tup in degrees_:
#     otxt += tup[0] + '\t' + str(tup[1]) + '\n'
    
# with open(f"degrees_{suffix}.txt", "w") as f:
#     f.write(otxt)




