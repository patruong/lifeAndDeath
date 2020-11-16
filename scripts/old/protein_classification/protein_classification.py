#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 14:45:41 2020

@author: ptruong
"""

import pandas as pd 
import numpy as np 

# Load files
df = pd.read_csv("peptide_tryptic_melted.csv", sep ="\t")
df = pd.read_csv("melted_treshold.csv", sep = "\t")

uniprot = pd.read_csv("uniprot-proteome UP000005640.tab", sep = "\t")
peptides = pd.read_csv("peptides.txt", sep = "\t")
peptides_tryptic = pd.read_csv("peptides_tryptic.txt", sep = "\t")


# Test to find GO from uniprot dB
protein = df.Leading_razor_protein[0]

def get_GO_from_uniprot_protein(protein):
    """
    protein - uniprot protein. E.g. "P55011"
    #uniprot - uniprot tab-seperated database with Gene ontology (GO) column.

    returns {"GO": <list of GO>, "number_of_GO":<int>}
    """
    res = uniprot[uniprot.Entry == protein]
    go = res["Gene ontology (GO)"].values
    if len(go) > 1:
        print("WARNING GO-value-array > 1")
    try:
        go = go[0].split(";")
    except:
        print(protein)
        result = {"GO":np.nan,"number_of_GO":np.nan}
        return 
    n_go = len(go)
    result = {"GO":go,"number_of_GO":n_go}
    return result

res = get_GO_from_uniprot_protein(protein)

df_col = df.Leading_razor_protein.apply(get_GO_from_uniprot_protein)

def adj_protein(protein):
    """
    Function maps proteins by removing suffix and prefix.
    
    Example of mapping.
    Q86U42 --> Q86U42
    Q16891-2 --> Q16891
    CON__P34955 --> P34955
    """
    split_1 = protein.split("-") # case with -
    split_2 = protein.split("__") # case with __
    if len(split_1) == 2: # case for protein with -
        protein = split_1[0]
        split_2 = protein.split("__")
        if len(split_2) == 2:
            protein = split_2[1]                
    elif len(split_2) == 2:
        protein = split_2[1]
    return protein

df["adj_protein"] = df.Leading_razor_protein.apply(adj_protein)

df_col = df.adj_protein.apply(get_GO_from_uniprot_protein)

adj_protein("CON__P01045-1")


 Q9NYF9

Type
Interactor (Protein)
Species
Homo sapiens
Synonyms
GOGA2_HUMAN, Q6GRM9, Q9BRB0, Q9NYF9, A0A0C4DGS5










