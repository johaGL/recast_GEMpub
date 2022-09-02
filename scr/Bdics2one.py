#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 18:58:41 2022

@author: johanna
examples:
python3 Bdics2one.py --mywdir ~/recast_GEMpub/ \
       --suffix brain_normal \
       --biomartfile "~/cocultureProj/biomart_db/bm_human.rds" 
       
python3 Bdics2one.py --mywdir ~/recast_GEMpub/ \
       --suffix GBM_tumor \
       --biomartfile "~/cocultureProj/biomart_db/bm_human.rds" 
   
NOTE: if in bb8 : --xml /scratch/CBIB/jgalvis/datapub/GEM/U-251MG.xml
"""
import os
import argparse
from fun_pickle2newpickle import *
from fun_bipartiteGEM import *

parser = argparse.ArgumentParser()
parser.add_argument("--xml" , help="full path to the xml")
parser.add_argument("--mywdir")
parser.add_argument("--suffix")
parser.add_argument("--biomartfile")
args = parser.parse_args()
print(args.xml)
print(args.suffix)


os.chdir(args.mywdir)
if not os.path.exists(os.path.join(os.getcwd(), args.suffix)):
    os.makedirs(args.suffix)  

"""    
B) use fun_pickle2newpickle
"""

myBM = pyreadr.read_r(args.biomartfile)

with open(f'{args.suffix}/reactions_{args.suffix}.pkl', 'rb') as f:
    reacsD = pk.load(f)   
    
with open(f'{args.suffix}/genes_metblits_{args.suffix}.pkl', 'rb') as f:
    genesmetblitsD = pk.load(f)
    
genesD, metsD = splitgemetdico(genesmetblitsD)
    
syD = ensemblsymD(myBM[None])  
    
outdf = genesdico2df(syD, genesD)
outdf.to_csv(f'{args.suffix}/gensymbols_{args.suffix}.csv')  

abool = isuniqval_inD(metsD, 'name')
print(f"    metabolites dictionary: {str(abool)}\n")
    # print(f"checking uniqueness of {akey} in dictionary")
    
print(showmultiplicitycause(metsD))

newmetsD = meltkeydico(metsD, "name", "compartment", "id")

enzy_reac_prod = doenzy_reac_prodII(reacsD, genesD, newmetsD, syD) 

print(f"   genes dictionary: {isuniqval_inD(genesD, 'name')}\n")
genesD = donewgenesD( outdf, genesD )
print(f"   genes new dictionary: {isuniqval_inD(genesD, 'symbol')}\n")
newgenes_sym_D = exchangekeydico(genesD, "symbol", "id")
#print(newgenes_sym_D["DPM1"])

# save result
groupedD = {'enzy_reac_prod':enzy_reac_prod,
            'genes' : newgenes_sym_D,
            'metabolites': newmetsD}  
with open(f'{args.suffix}/dictioFull_{args.suffix}.pkl','wb') as f:
    pk.dump(obj=groupedD, file=f)

print("saved dictionary of dictionaries: enzy_reac_prod, genes and metabolites")
print(f'\t {args.suffix}/dictioFull_{args.suffix}.pkl')
        
"""
C)  bipartite
"""    
with open(f'{args.suffix}/dictioFull_{args.suffix}.pkl', "rb") as f:
    GEM = pickle.load(f)  
print(GEM.keys())
edges_ = allmetaboedges(GEM)

edges4txt = [ "\t".join(i)+"\n" for i in edges_ ]

bipafile = f'{args.suffix}/bipartite_{args.suffix}.txt'
with open(bipafile, "w") as f:
    f.writelines(edges4txt)
    
sms = doreportbipartite(args.suffix,GEM, edges_)
with open(f"{args.suffix}/reportbipartite_{args.suffix}.txt", "w") as f:
    f.writelines(sms)    
    
G = nx.Graph()
G.add_edges_from(edges_)
boolbipartite = nx.is_bipartite(G)
print(f"The graph from these Edges, is it bipartite?: *{boolbipartite}*")
# add line to auto report:
with open(f"{args.suffix}/reportbipartite_{args.suffix}.txt", "a") as f:
    f.writelines([f"\nis bipartite : {boolbipartite}"])

print(f"saved bipartite graph: {bipafile} ")
