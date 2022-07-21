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
df = pd.DataFrame(degrees_, columns= ['vertex', 'degree'])

for i, ro in df.iterrows():
    if "utamate|" in ro['vertex']:
        print(ro['degree'])

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
ax1.set_xlabel("....")
ax1.set_title(f'')
plt.savefig("mydegrees.txt")

init_q = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,  0.9, 0.99,1.0}
quantiles = df['degree'].quantile(list(init_q))
quantiles

above_09 = df.loc[df['degree'] >= quantiles[0.99]]

# 
# otxt = ''
# for tup in degrees_:
#     otxt += tup[0] + '\t' + str(tup[1]) + '\n'
    
# with open(f"degrees_{suffix}.txt", "w") as f:
#     f.write(otxt)




