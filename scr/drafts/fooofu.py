#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ugly quasi-copy  of GEMformatter 29 april, ok tests

@author: johanna
examples:
python3 ????.py --xml ~/datapub/GEM/0_dwld/brain.xml \
       --mywdir ~/recast_GEMpub/ --suffix brain_normal \
       --biomartfile "~/cocultureProj/biomart_db/bm_human.rds" 
       
python3 ????.py --xml ~/datapub/GEM/0_dwld/U-251MG.xml \
   --mywdir ~/recast_GEMpub/ --suffix GBM_tumor \
   --biomartfile "~/cocultureProj/biomart_db/bm_human.rds" 
   

"""
import os
import sys
import argparse
from fun_GEMxml2pickle import *
from fun_pickle2newpickle import *

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

# # use GEMxml2pickle functions
# thexml = args.xml
# docnodes, docedges = xml_doc_parsing(args.xml)
# readico = set_reaction_dict(docedges) 
# complexdico = set_metabo_dict(docnodes)

# with open(f'{args.suffix}/reactions_{args.suffix}.pkl', 'wb') as f:
#     pk.dump( readico, f )
# with open(f'{args.suffix}/genes_metblits_{args.suffix}.pkl', "wb") as f:
#     pk.dump( complexdico, f )
     
#use2_pickle2newpickle

myBM = pyreadr.read_r(args.biomartfile)

with open(f'{args.suffix}/reactions_{args.suffix}.pkl', 'rb') as f:
    reacsD = pk.load(f)   
    
with open(f'{args.suffix}/genes_metblits_{args.suffix}.pkl', 'rb') as f:
    genesmetblitsD = pk.load(f)
    
genesD, metsD = splitgemetdico(genesmetblitsD)
print("..........");print(genesD['E_2570'])    
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
newgenesD = donewgenesD( outdf, genesD )
print(f"   genes new dictionary: {isuniqval_inD(newgenesD, 'symbol')}\n")
print(genesD['E_2570'])
print(genesD['ENSG00000168389'])
newgenes_sym_D = exchangekeydico(newgenesD, "symbol", "id")
#print(newgenes_sym_D["DPM1"])


# save result
groupedD = {'enzy_reac_prod':enzy_reac_prod,
            'genes' : newgenes_sym_D,
            'metabolites': newmetsD}  
with open(f'{args.suffix}/HMAdictio_{args.suffix}.pkl','wb') as f:
    pk.dump(obj=groupedD, file=f)

print("saved dictionary of dictionaries: reactions, enzymes and metabolites")
print(f'\t {args.suffix}/HMAdictio_{args.suffix}.pkl')
        
     