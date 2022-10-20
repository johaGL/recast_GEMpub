#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:01:51 2022

From newest json human GEM:
    which was downloaded from:
        https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.yml
    note: SysBioChalmers is the team of the Human Metabolic Atlas

obtain two types of edges files:
    Bipartite
    and 
    not bipartite with edges:
        - joining metabolites
        - having a list of enzymes as decorator
@author: johanna

usage
python3 2022_gem_human.py \
    --yml ~/recast_GEMpub/data/GEM2022/Human-GEM.yml \
    --mywdir ~/recast_GEMpub/ --suffix hum_2022 \
        --biomartfile "~/cocultureProj/biomart_db/bm_human.rds"  
        
   
"""


import os
import argparse
from fun_pickle2newpickle import *
from fun_bipartiteGEM import *
import yaml
import re
import pandas as pd

def doc2rsrs(doc):
    """
    create dico REACTION_ID : {'enzymes': [xx|c, yy|...], 'metabolites' : []}  
    """
    rsrs = dict()
    allgenesgem = set()        
    for k_reaction in range(len(doc[2][1])):  # len(doc[2][1])
        llhere = doc[2][1][k_reaction]       
        k_genes = re.split(" or | and ", llhere[5][1])
        k_genes = [i.replace("(", "").replace(")","") for i in k_genes]
        if k_genes != ['']: # no not add if not gene associated
            #print(llhere)
            #print(k_genes)        
            allgenesgem.update(k_genes)            
            #print(llhere)
            k_id_react = llhere[0][1]
            k_mets = [i for (i,j) in llhere[2][1]]
            #print(k_mets)            
            if k_id_react not in rsrs.keys():
                rsrs[k_id_react] = { 'enzymes' : k_genes,
                                    'metabolites' : k_mets}
            else:
                print(f"ERROR, this reaction {k_id_react} is repeated")
    return rsrs, allgenesgem


def replaceidsbysymbol(namesbefore, syD):
    """
    Replace ensembl by symbols using biomart from my local:
    Parameters
    ----------
    namesbefore : list
    syD : dict

    Returns
    -------
    namesafter : list
    """
    namesafter = []
    for nami in namesbefore:
        try:
            namesafter.append(syD[nami])
        except:
            namesafter.append(nami) # ENSG if failed find symbol
    return namesafter


def update_rsrs_genes(Xbm, allgenesgem, rsrs):        
    # check matches : 
        # use these two functions from fun_pickle2newpickle
    myBM = pyreadr.read_r(os.path.expanduser(Xbm))
    syD = ensemblsymD(myBM[None])
    nono = list()
    for i in list(allgenesgem):
        try: 
            fii = syD[i]
            if fii.startswith("ENSG"):
                print(fii)
        except:
            nono.append(i)
    print("not found : ",  nono, " \n>>left as it<<")
    # do replacement 
    for k_reac in rsrs.keys():
        genesbefore = rsrs[k_reac]['enzymes']
        genesafter = replaceidsbysymbol(genesbefore, syD)
        rsrs[k_reac]['enzymes'] = genesafter
    return rsrs
 
    
def dometsdico(doc):
    """
    search metabolites names in yaml itself, do dictio
    new ! : only compartment e (extracellular) kept,
            all the other replaced by c (intracellular)
    """
    mets2022 = dict()
    k_met = 0
    for k_met in range(len(doc[1][1])):        
        doc[1][1][k_met]
        # [('id', 'MAM00001c'),
        # ('name', '(-)-trans-carveol'),
        # ('compartment', 'c'),
        # ('formula', 'C10H16O'), ... ]
        k_id_met = doc[1][1][k_met][0][1]
        k_name_met = doc[1][1][k_met][1][1] 
        k_compartment_met = doc[1][1][k_met][2][1]
        if k_id_met not in mets2022.keys():
            if k_compartment_met == "e": # new !! 
                mets2022[k_id_met] = k_name_met +"|"+k_compartment_met
            else:
                mets2022[k_id_met] = k_name_met +"|ic" 
        else:
            "ERRORRRRRR: repeated metabolite id" 
    return mets2022    
 
   
def removepromiscuous(metabolites_, promiscuous_):
    elementstodrop = list()
    for i in metabolites_:
        if i.split("|")[0] in promiscuous_:
            elementstodrop.append(i)
    outm = list(set(metabolites_) - set(elementstodrop))
    return outm
    
def updatemets(rsrs, mets2022):    
    # do replacement: the metabolites names
    for k_reac in rsrs.keys():
        metsbefore = rsrs[k_reac]['metabolites']
        metsafter = replaceidsbysymbol(metsbefore, mets2022)
        
        rsrs[k_reac]['metabolites'] = metsafter
    return rsrs



def rsrs2bipatxt(rsrs, ofile, promiscuous):
    alleds = list()
    for k_reac in rsrs.keys():
        for c in rsrs[k_reac]['enzymes']:
            metsclean = removepromiscuous(rsrs[k_reac]['metabolites'], promiscuous)
            for d in metsclean:
                tmp = c +"\t" + d + "\n"
                if tmp not in alleds:
                    alleds.append(tmp)
    with open(ofile, "w") as f:
        f.writelines(alleds)
    return 0

def rsrs2metedges(rsrs, ofile, promiscuous):
    litups = dict()
    for k_reac in rsrs.keys():
        metsclean = removepromiscuous(rsrs[k_reac]['metabolites'], promiscuous)
        slici = 0
        while slici < len(metsclean):
            x = metsclean[slici]
            for y in metsclean[(slici+1): ]:
                if x != y :                   
                    if tuple(sorted([x,y])) not in litups.keys():
                        litups[tuple(sorted([x,y]))] = rsrs[k_reac]["enzymes"]
                    else:
                        litups[tuple(sorted([x,y]))] += rsrs[k_reac]["enzymes"]
            slici += 1
    alleds = list()
    for ktup in litups.keys():
        print(litups[ktup])
        enzymesjoined = " ".join(litups[ktup])  
        x = ktup[0]
        y  = ktup[1]
        tmp = f"{x}\t{y}\t{enzymesjoined}\n" 
        alleds.append(tmp)
    with open(ofile, "w") as f:
        f.writelines(alleds)
    return 0
            

#############
###### START
#############
parser = argparse.ArgumentParser()
parser.add_argument("--yml" , help="full path to the YML")
parser.add_argument("--mywdir")
parser.add_argument("--suffix")
parser.add_argument("--biomartfile")
args = parser.parse_args()
print(args.suffix)


promiscuous = ["H2O", "H+", "H2O2", "O2", "Na+", "CO2" , "CO",
               "[protein]", "phosphoprotein", "Pi", "PPi",
               "AMP", "ADP" , "ATP",
               "CMP", "CDP",  "GMP", "GDP", "UMP", "UDP",
               "NAD+","NADH", "NADP+", "NADPH", 
                "CoA", "acetyl-CoA" ] 

os.chdir(args.mywdir)
if not os.path.exists(os.path.join(os.getcwd(), args.suffix)):
    os.makedirs(args.suffix)  
    

print("Load yaml")
yf = os.path.expanduser(args.yml)
with open(yf, 'r') as f:
    doc = yaml.load(f, Loader=yaml.Loader)
print("parsing and building edges")
 #### use functions   
rsrs, allgenesgem = doc2rsrs(doc)
len(rsrs)
len(allgenesgem)  
print("end step 1 ")
rsrs = update_rsrs_genes(args.biomartfile, allgenesgem, rsrs)
# METABOLITES part                
mets2022  = dometsdico(doc)  
rsrs = updatemets(rsrs,mets2022)
print("end step 2 ")

ofilebip = f"{args.suffix}/bipatrimmed_{args.suffix}.txt"
rsrs2bipatxt(rsrs, ofilebip, promiscuous)
ofilemets = f"{args.suffix}/metsmets_{args.suffix}.txt"
rsrs2metedges(rsrs, ofilemets, promiscuous)


# print compartments new nomenclature 2022
cdf = pd.DataFrame(doc[4][1], columns = ["abbreviation", "compartment"])
cdf.to_csv(f"{args.suffix}/compartmentslegend_{args.suffix}.txt", index = False) 