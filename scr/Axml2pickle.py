#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 18:58:41 2022

@author: johanna
examples:
python3 Axml2pickle.py --xml ~/datapub/GEM/brain.xml \
       --mywdir ~/recast_GEMpub/ --suffix brain_normal
       
python3 Axml2pickle.py --xml ~/datapub/GEM/U-251MG.xml \
   --mywdir ~/recast_GEMpub/ --suffix GBM_tumor
   
NOTE: if in bb8 : --xml /scratch/CBIB/jgalvis/datapub/GEM/U-251MG.xml
"""
import os
import argparse
from fun_GEMxml2pickle import *

parser = argparse.ArgumentParser()
parser.add_argument("--xml" , help="full path to the xml")
parser.add_argument("--mywdir")
parser.add_argument("--suffix")
args = parser.parse_args()
print(args.xml)
print(args.suffix)


os.chdir(args.mywdir)
if not os.path.exists(os.path.join(os.getcwd(), args.suffix)):
    os.makedirs(args.suffix)  
"""
A) fun_GEMxml2pickle functions (commented because is loong, already done)
"""

print(f"parsing xml reactions and metabolites ")

docnodes, docedges = xml_doc_parsing(args.xml)
readico = set_reaction_dict(docedges) 
complexdico = set_metabo_dict(docnodes)

print(f"saving reactions and metabolites to pickle files at : \
      {args.mywdir} ==> {args.suffix}/")
with open(f'{args.suffix}/reactions_{args.suffix}.pkl', 'wb') as f:
    pk.dump( readico, f )
with open(f'{args.suffix}/genes_metblits_{args.suffix}.pkl', "wb") as f:
    pk.dump( complexdico, f )

print("DONE")