#!/usr/bin/env python3
"""
compare edges of these two different sets:
union brain tumor
vs
total human GEM 2022
"""

import os
import pandas as pd

def opentxtedges(filetxt):
    """
    filetxt : a tabular delimited file of edges
    lotu_ = list of tuples
    """
    with open(filetxt, 'r') as f:
        tmp = f.readlines()
    lotu_ = [tuple(i.strip().split('\t')) for i in tmp]
    return lotu_

def detectdifference_transporters(onetup_, twotup_):
    # SLC16 and SLC1A families to detect:
    def countedgesmatching_solutecarrier(anytup_):
        tmp = pd.DataFrame(anytup_)
        tmp = tmp.dropna()
        k = tmp.loc[tmp[0].str.contains('SLC16')]
        l = tmp.loc[tmp[1].str.contains('SLC16')]
        m = tmp.loc[tmp[0].str.contains('SLC1A')]
        n = tmp.loc[tmp[1].str.contains('SLC1A')]
        x = k.shape[0] + l.shape[0] + m.shape[0] + n.shape[0]
        return x
    oneq = countedgesmatching_solutecarrier(onetup_)
    twoq = countedgesmatching_solutecarrier(twotup_)
    return oneq, twoq

def suppresscompartment(lotu_):
    newout = list()
    for tup in lotu_:
        # as any of both nodes in tup can be the metabolite:
        if len(tup) == 2 and (tup[0] != '') and (tup[1] != ''):
            i = tup[0].split("|")[0]
            j = tup[1].split("|")[0]
        # hey, discovered this error : edges with one single element:
        if len(tup) < 2:
            # print(tup) # pending to delete this type of singlets
            continue
    newout.append((i,j))
    return newout

humtotal = "hum_whole/hum_whole_trimmed.txt"
union_brain_GBM = "union_brain_GBM/union_brain_normal_GBM_tumor.txt"
myw = os.path.expanduser("~/recast_GEMpub/")

humtotu_ = opentxtedges(myw + humtotal)
uniotu_ = opentxtedges(myw + union_brain_GBM)

humq, unioq = detectdifference_transporters(humtotu_, uniotu_)

print("\nComparing : union brain tumor -vs- total human GEM 2022 ")

print("\nnumber of edges containing lactate or glutamate transporters" )
print("in total human GEM :", humq)
print("in union brain tumor:", unioq)

print("\nimportant! : as compartment nomenclature is different (it changed in 2022)\n \
  just use metabolites without their compartment (number of edges preserved in both sides)")

humtotu_ = suppresscompartment(humtotu_)
uniotu_ = suppresscompartment(uniotu_)

print("\nDifferences among versions, this time with suppressed compartment:\n")
print(" ->  present in human-gem, but not in union-brain-tumor:")
print(set(humtotu_) - set(uniotu_)) # tup[0].split("|")[0]
print(" ->  present in union-brain-tumor, but not in human-gem:")
print(set(uniotu_) - set(humtotu_)) # {('ALDH9A1', '4-aminobutyrate')}

m = "\nin conclusion, differences are negligible after suppressing the compartments"
print(m)