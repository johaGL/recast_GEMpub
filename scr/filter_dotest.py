#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 17:47:26 2022

@author: johanna
"""

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# parser = argparse.ArgumentParser()
# parser.add_argument("--mywdir")


#os.chdir(args.mywdir)

mywdir = "~/recast_GEMpub/"
os.chdir(os.path.expanduser(mywdir))

bipafile = "brain_normal/bipartite_brain_normal.txt"
suffix = "brain_normal"

dfi = pd.read_csv(bipafile, header=None, sep="\t")

G = nx.Graph()

for i, ro in dfi.iterrows():
    G.add_edge(ro[0], ro[1])
    
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
ax1.set_title(f'')
plt.savefig("mydegrees.pdf")

init_q = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.99, 1.0}
quantiles = df['degree'].quantile(list(init_q))
quantiles
print("\n OPTION 1 : Using a percentile threshold")

chosenpercentile = 0.99
above_perc = df.loc[df['degree'] >= quantiles[chosenpercentile]]
print(f"printing vertices whose degree is found above percentile {chosenpercentile} ")
print(above_perc)
print("Note that some interesting genes are in risk of being excluded if using percentile cutoff")

# another method:
unaccep_degree = 0
print(f"\n OPTION 2 : \n\
      a) Using a list of promiscuous vertices, and then \n\
      b) remove  vertices with degree <= {unaccep_degree}\n")

promiscuous = ["H2O", "H+", "H2O2", "O2", "Na+", "CO2" , "CO",
               "[protein]", "phosphoprotein", "Pi", "PPi",
               "AMP", "ADP" , "ATP",
               "CMP", "CDP",  "GMP", "GDP", "UMP", "UDP",
               "NAD+","NADH", "NADP+", "NADPH", 
                "CoA", "acetyl-CoA" ] 

print(f"Deleting promiscuous vertices {promiscuous} \n \
      (all compartments comprised), and all adjacent edges.\n")
      

Gm = G.copy()
preexistingnodes = list(Gm.nodes)
promiscuous_wicompartment = set([ n for n in preexistingnodes if n.split("|")[0] in promiscuous])
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
    print(condicount)
    print(df_m_nogood)

print(f'also nodes with degree <= {unaccep_degree} removed:\n {suppressedall}')
for i in suppressedall:
    if i.startswith("SLC16"):
        print(i)
        
len(suppressedall)        
        
degrees_m_ = list(Gm.degree())
df_m = pd.DataFrame(degrees_m_, columns = ['vertex', 'degree'])


print("verifying new degrees")
print(df_m.sort_values("degree", axis=0, ascending=False))


print(df_m.shape)



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




