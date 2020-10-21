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
