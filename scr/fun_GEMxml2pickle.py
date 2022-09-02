#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 19:46:23 2022

@authors: A Coffe, johaGL
"""

from xml.dom import minidom as mdom
import pickle as pk
import os

def set_metabo_dict(xml_nodes):
    metabo_dict = {}
    for s in xml_nodes:
        attr = s.attributes
        species = dict(attr.items())
        n = species.pop('id')
        metabo_dict[n] = species

    return metabo_dict


def extract_reactant_products(xml_listof):
    mol_stack = []
    if xml_listof is not None:
        for i, mol in enumerate(xml_listof):
            tmp = dict(mol.attributes.items())
            if len(tmp) == 0:
                for child in mol.childNodes:
                    tmp = dict(child.attributes.items())
                    mol_stack.append(tmp.get('species'))
            else:
                mol_stack.append(tmp.get('species'))

    return mol_stack


def xml_doc_parsing(path):
    """

    :param path: xml species/reaction  list pathway
    :return:
    """
    myxml = mdom.parse(path)
    xml_nodes = myxml.getElementsByTagName('species')
    xml_edges = myxml.getElementsByTagName('reaction')

    return xml_nodes, xml_edges

def set_reaction_dict(xml_edges):
    reaction_dict = {}
    for r in xml_edges:
        ids = dict(r.attributes.items())
        notes = r.getElementsByTagName('notes')             
       
        listOfReactants = r.getElementsByTagName('listOfReactants')[0].childNodes       
        listOfProducts = r.getElementsByTagName('listOfProducts')     
        listOfModifiers = r.getElementsByTagName('listOfModifiers')
        n = ids.pop('id')
        ids.pop('name')

        reactants = extract_reactant_products(listOfReactants)
        
        products = extract_reactant_products(listOfProducts)
        supp = notes[0].lastChild.childNodes[0].childNodes[0].data
        E_ =  extract_reactant_products(listOfModifiers)
       
        edge_list = []

        for o in reactants:  # this was proposed by A coffe
            for t in products:
                w = (o, t)
                edge_list.append(w)
            
            dkey = n + '_' + ids.get('metaid')
            reaction_dict[dkey] = {'id': ids, 
                                   'modifiers' : E_,
                                   'reaction': edge_list, 
                                   'reactants' : reactants,
                                    'products' : products,
                                   'notes': supp}

    return reaction_dict

if __name__ == "__main__":
    thexml = os.path.expanduser("~/datapub/GEM/0_dwld/brain.xml") 
    suffix = "brain_normal"
    mywdir = os.path.expanduser("~/recast_GEMpub/")
    os.chdir(mywdir)
    os.makedirs(suffix)
    
    suffix = "brain_normal"
    print(thexml)
    docnodes, docedges = xml_doc_parsing(thexml)
    readico = set_reaction_dict(docedges) 
    complexdico = set_metabo_dict(docnodes)
    
    with open(f'{suffix}/reactions_{suffix}.pkl', 'wb') as f:
        pk.dump( readico, f )
    with open(f'{suffix}/genes_metblits_{suffix}.pkl', "wb") as f:
        pk.dump( complexdico, f )
