#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 11:31:12 2022

Generates pickle dictionaries of the newest full human GEM (21 june 2022)
which was downloaded from:
    https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.yml
note: SysBioChalmers is the team of the Human Metabolic Atlas
     
NEW observed features compared to the ancient GEMs : 
    - no "modifiers" "reactants" "products" keywords
    - compartments are VERY different, they are printed with my code here
    
Usage:
  
python3 newHumanGEMformatter.py \
    --yml ~/recast_GEMpub/data/GEM2022/Human-GEM.yml \
    --mywdir ~/recast_GEMpub/ --suffix hum_whole \
        --biomartfile "~/cocultureProj/biomart_db/bm_human.rds"    
        
        
 @author: johanna  
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
            mets2022[k_id_met] = k_name_met +"|"+k_compartment_met
        else:
            "ERRORRRRRR: repeated metabolite id" 
    return mets2022    
 
    
def updatemets(rsrs, mets2022):    
    # do replacement: the metabolites names
    for k_reac in rsrs.keys():
        metsbefore = rsrs[k_reac]['metabolites']
        metsafter = replaceidsbysymbol(metsbefore, mets2022)
        rsrs[k_reac]['metabolites'] = metsafter
    return rsrs


def rsrs2bipatxt(rsrs, ofile):
    otxt = ''
    for k_reac in rsrs.keys():
        for c in rsrs[k_reac]['enzymes']:
            for d in rsrs[k_reac]['metabolites']:
                tmp = c +"\t" + d + "\n"
                otxt += tmp
    with open(ofile, "w") as f:
        f.write(otxt)
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

ofile = f"{args.suffix}/{args.suffix}_bipartite.txt"
rsrs2bipatxt(rsrs, ofile)


# print compartments new nomenclature 2022
cdf = pd.DataFrame(doc[4][1], columns = ["abbreviation", "compartment"])
cdf.to_csv(f"{args.suffix}/compartmentslegend_{args.suffix}.txt", index = False) 

# global report
with open(ofile, "r") as f:
    txtedg = f.readlines()

edges_ = [tuple(i.strip().split("\t")) for i in txtedg]
#sms = doreportbipartite(args.suffix, rsrs, edges_) # failed, keys 2022 different
A = f"---------------------------------------- \
       \n{args.suffix}\n----------------------------------------"
a = "\nTransformed GEM dictionnary into list of Edges\n"
b = f"number of reactions *involving gene(s)* : {len(rsrs)}\n"
c = f"Edges : {len(edges_)}\
  \nNodes\
  \n\tgenes : {len(allgenesgem)}\n\tmetabolites : {len(mets2022)}"
sms = [A, a, b, c]
with open(f"{args.suffix}/reportbipartite_{args.suffix}.txt", "w") as f:
    f.writelines(sms)
    
print("END")

# END exec
 
# if __name__ == "__main__":
#     Xbm = "~/cocultureProj/biomart_db/bm_human.rds" 
#     yf = os.path.expanduser("~/Human-GEM.yml")
#     yf
#     suffix = "hum_whole"
#     with open(yf, 'r') as f:
#         doc = yaml.load(f, Loader=yaml.Loader)
        
#     for i in range(len(doc)):
#         print(i, doc[i][0])
#         print("     ", len(doc[i][1]))
#         #print(i, len(doc[i][2]))
#     # 0 metaData
#     #       10
#     # 1 metabolites
#     #       8369
#     # 2 reactions
#     #       13070
#     # 3 genes
#     #       3067
#     # 4 compartments
#     #       9  
   
#     # One single reaction view:
#     for k in range(len(doc[2][1][0])):
#         print(k, doc[2][1][0][k])
#     # 0 ('id', 'MAR03905')
#     # 1 ('name', 'ethanol:NAD+ oxidoreductase')
#     # 2 ('metabolites', [('MAM01249c', 1), ('MAM01796c', -1), ('MAM02039c', 1), ('MAM02552c', -1), ('MAM02553c', 1)])
#     # 3 ('lower_bound', 0)
#     # 4 ('upper_bound', 1000)
#     # 5 ('gene_reaction_rule', 'ENSG00000147576 or ENSG00000172955 or ENSG00000180011 or ENSG00000187758 or ENSG00000196344 or ENSG00000196616 or ENSG00000197894 or ENSG00000198099 or ENSG00000248144')
#     # 6 ('rxnNotes', '')
#     # 7 ('rxnFrom', 'HMRdatabase')
#     # 8 ('eccodes', '1.1.1.1;1.1.1.71')
#     # 9 ('references', 'PMID:10868354;PMID:12491384;PMID:12818203;PMID:14674758;PMID:15289102;PMID:15299346;PMID:15327949;PMID:15682493;PMID:15713978')
#     # 10 ('subsystem', ['Glycolysis / Gluconeogenesis'])
#     # 11 ('confidence_score', 0)
  
#      #### use functions   
#     rsrs, allgenesgem = doc2rsrs(doc)
#     len(rsrs)
#     len(allgenesgem)  
#     rsrs = update_rsrs_genes(Xbm, allgenesgem, rsrs)
#     # METABOLITES part                
#     mets2022  = dometsdico(doc)  
#     rsrs = updatemets(rsrs,mets2022)
    
 

    
