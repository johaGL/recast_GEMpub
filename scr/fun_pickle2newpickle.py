#!/usr/bin/env python3
"""
Created on Thu Apr  7 14:32:58 2022
@author: johaGL

Needs  :  EXAMPLE
    suffix : "brain_normal"
    localpathBM : full path to local biomart dataframe as .rds file
    awdir : "~/recast_GEMpub/"
"""

import pickle as pk
import os
import pandas as pd
import pyreadr

def splitgemetdico(genesmetblitsD):
    """
    split genes_metblits dictionary into 2 dictionaries:
        - genes and 
        - metabolites/cofactors 
    """
    genesD = {}
    metsD = {}
    for uk in genesmetblitsD.keys():
        if genesmetblitsD[uk]['name'].startswith("ENSG"):
            genesD[uk] = genesmetblitsD[uk] # genes
        else: 
            metsD[uk] = genesmetblitsD[uk] # metaboltes
    print("end splitgemetdico")
    return genesD, metsD

def ensemblsymD(BMdf):
    """
    get symbols for each ensembl id in BMdf dataframe
    BMdf dataframe has columns "ensembl_gene_id" and "external_gene_name"
          as given by biomart
    """
    assert isinstance(BMdf, pd.core.frame.DataFrame), \
        "ERROR! must be pandas dataframe input for ensemblsymD(BMdf)"
    syD = {}
    for i, row in BMdf.iterrows():
        ensid = str(row["ensembl_gene_id"])
        symbo = str(row["external_gene_name"])
        syD[ensid] = symbo
    print("finished dico from BioMart df/n")
    #print(BMdf[BMdf['external_gene_name']== "PHOSPHO2"]["ensembl_gene_id"])
    return syD
    

def genesdico2df(syD, genesD):
    """
    get ensembl and symbols dataframe:
    """  
    ensids = []
    symbols = []; notfoundcount=0;
    for k in genesD.keys():
        ensids.append(genesD[k]['name'])
        try:
            symbols.append(syD[ genesD[k]['name'] ])
        except:
            print(f"{genesD[k]['name']} ===> NOT FOUND")
            notfoundcount += 1
            symbols.append(genesD[k]['name'])
            
    outdf = pd.DataFrame(data={'ensid': ensids,
                               'symbol': symbols})
    print(f'successfully found {len(symbols) - notfoundcount} gene symbols.\
          \nThe NOTFOUND ({notfoundcount}) were assigned ENSG000...\n')
    return outdf

def doenzy_reac_prodII(reacsD, genesD, newmetsD, syD):
    """
    version II: do enzy_reac_prod dictionary: reactants, gene symbols, products
    """
    tmpmeD = exchangekeydico(newmetsD, 'id', 'name_full') 
    # exchangekeydico: "urate|C_x" to 'name_full' key,  id "M_m03120x" main key
    enzy_reac_prod = {}
    count = 0
    for k in reacsD.keys():
        reactants_ = []
        products_ = []
        enzymes_ = []
        for h in reacsD[k]['reactants']:
            reactants_.append( tmpmeD[ h ]['name_full'] )
        for h in reacsD[k]['products']:
            products_.append( tmpmeD[h]['name_full'] )
        for h in reacsD[k]['modifiers']:
            ensid = genesD[h]['name']
            try:
                mysymbol = syD[ ensid ]
                enzymes_.append(mysymbol)
            except: 
                continue
        if len(enzymes_) > 0 :    
            enzy_reac_prod[count] ={ 'enzymes' : enzymes_,
                                    'reactants' : reactants_,
                                    'products' : products_
                                    }
            count += 1 
    return enzy_reac_prod

def doenzy_reac_prod(reacsD, genesD, metsD, symD):
    """
    do enzy_reac_prod dictionary: reactants, gene symbols, products
    """
    enzy_reac_prod = {}
    count = 0
    for k in reacsD.keys():
        reactants_ = []
        products_ = []
        enzymes_ = []
        for h in reacsD[k]['reactants']:
            reactants_.append( metsD[ h ]['name'] )
        for h in reacsD[k]['products']:
            products_.append( metsD[h]['name'] )
        for h in reacsD[k]['modifiers']:
            ensid = genesD[h]['name']
            try:
                mysymbol = symD[ ensid ]
                enzymes_.append(mysymbol)
            except: 
                continue
        if len(enzymes_) > 0 :    
            enzy_reac_prod[count] ={ 'enzymes' : enzymes_,
                                    'reactants' : reactants_,
                                    'products' : products_
                                    }
            count += 1 
    return enzy_reac_prod 

def isuniqval_inD(nestedD, akey):
    """
    Parameters
    ----------
    nestedD :  dictionary
        nested dictionary, two or more levels
    akey : str
        second level key, ex: d[k]["name"] , akey = "name"

    Returns
    -------
    boolean
    """  
    print(f"checking uniqueness of {akey} in dictionary")
    allvals = []
    for k in nestedD.keys():
        allvals.append(nestedD[k][akey])
    if (len(allvals) == len(set(allvals))):
        return True
    else: 
        return False
    

def donewgenesD(outdf, genesD):   
    """
    Parameters
    ----------
    outdf : pandas df
    genesD : dictionary

    Returns
    -------
    newgenD : dictionary, with:
        key level1 "ENSG...", 
        sub-keys metaid, name, compartment, initialAmount, sboTerm
             C_c,     symbol,         id (id is "E_.." (suffix of metaid)) 
    """
    newgenD = {}
    for k in genesD.keys():        
        fee = genesD[k]['name']
        #symbo = outdf[ outdf['ensid'] == fee]["symbol"].to_string().split()[1]
        row = outdf[ outdf['ensid'] == fee]
        symbo = row.iloc[0]["symbol"]
        newgenD[ fee ] = genesD[k]
        newgenD[ fee ]['symbol'] = symbo
        newgenD[ fee ]['id'] = k
    print("finished new genes dictionary, with symbol subkey added")
    return newgenD

def exchangekeydico(nestedD, newl1key,newl2key):
    """
    nestedD : dictionary 
    newl1key :level one key
    newl2key : level two key, must be unique as newl1key is
        DESCRIPTION.
    Returns : dictionary with keys exchanged
    """
    newD = {}
    if isuniqval_inD(nestedD, newl1key):
        for k in nestedD.keys():
            nestedD[k][newl2key] = k
            newD[nestedD[k][newl1key]] = nestedD[k]
        return newD
    else:
        print(f"impossible to set {newl1key} as newkey, it is not unique")
        return "ERROR newkey not unique"

def showmultiplicitycause(metsD):
    repeats = {}
    for k in metsD.keys():
        ametabli = metsD[k]["name"]
        if ametabli not in repeats.keys():
            repeats[ametabli] = 1
        else:
            repeats[ametabli] += 1
    print(" * each metaid is a unique association of id and compartment :*")
    print(f"    showing max case ===> {max(repeats, key=repeats.get )} <=== :")
    for k in metsD.keys():
        if metsD[k]["name"] == max(repeats, key=repeats.get):
            print(f"     *    {k}    : {metsD[k]}*")
    return 0

def readorpipe(mysep):
    """"returns a pipe if no separator character specified"""
    return mysep or "|"

def meltkeydico(nestedD, keyl2a, keyl2b, keytodown, mysep=None):
    SEP = readorpipe(mysep) # by default returns "|"
    newD = {}
    for k in nestedD.keys():
        meltedkey = f'{nestedD[k][keyl2a]}{SEP}{nestedD[k][keyl2b]}'
        newD[meltedkey] = nestedD[k]
        newD[meltedkey][keytodown] = k
    return newD
        
    
if __name__=='__main__':
    
    suffix = "GBM_tumor"
    localpathBM = "~/cocultureProj/biomart_db/bm_human.rds" # local biomart 
    mywdir = os.path.expanduser("~/recast_GEMpub/")
    
    myBM = pyreadr.read_r(os.path.expanduser(localpathBM))
    os.chdir(mywdir)

    with open(f'{suffix}/reactions_{suffix}.pkl', 'rb') as f:
        reacsD = pk.load(f)   
        
    with open(f'{suffix}/genes_metblits_{suffix}.pkl', 'rb') as f:
        genesmetblitsD = pk.load(f)
        
    genesD, metsD = splitgemetdico(genesmetblitsD)
        
    syD = ensemblsymD(myBM[None])  
        
    outdf = genesdico2df(syD, genesD)
    outdf.to_csv(f'{suffix}/gensymbols_{suffix}.csv')  
    
    print(f"   genes dictionary: {isuniqval_inD(genesD, 'name')}")
    newgenesD = donewgenesD( outdf, genesD )
    print(f"   genes new dictionary: {isuniqval_inD(newgenesD, 'symbol')}")
    
    newgenes_sym_D = exchangekeydico(genesD, "symbol", "id")
    print(newgenes_sym_D["DPM1"])
    print(f"    metabolites dictionary: {isuniqval_inD(metsD, 'name')}")
        # print(f"checking uniqueness of {akey} in dictionary")
        
    print(isuniqval_inD(metsD, "name"))
    print(showmultiplicitycause(metsD))
    newmetsD = meltkeydico(metsD, "name", "compartment", "id")
    
    enzy_reac_prod = doenzy_reac_prodII(reacsD, genesD, newmetsD, syD)    

    # save result
    groupedD = {'enzy_reac_prod':enzy_reac_prod,
                'genes' : newgenes_sym_D,
                'metabolites': newmetsD}  
    with open(f'{suffix}/dictioFull_{suffix}.pkl','wb') as f:
        pk.dump(obj=groupedD, file=f)
            
    ## END
    
    # crossrefs_ = [] 
    # systems_ = []
    # for rk in reacsD.keys():
    #     print(reacsD[rk])
    #     crossrefs_.append(reacsD[rk]['modifiers'])
    #     systems_.append(reacsD[rk]['notes'])
    # print(len(crossrefs_))
    #print(set(systems_))
    
    # ...
    
    # print(reacsD["R_HMR_1_metaid_R_HMR_1"])   
    # themets =  []
    # for tup in reacsD["R_HMR_1_metaid_R_HMR_1"]['reaction']:
    #     themets.append(tup[0])
    #     themets.append(tup[1])
    # print(set(themets))
    # print("---- ? ----")
    # for amet in themets:
    #     print(metsD[amet])
    # # the enzyme:
    # print(genesD["E_1"])
    
    
    
