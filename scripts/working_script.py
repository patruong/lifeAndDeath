#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 13:56:52 2021

@author: patrick
"""

import pandas as pd

peptides_tryptic = pd.read_csv("peptides tryptic.txt", sep = "\t")
protein_tryptic = pd.read_csv("proteinGroups tryptic.txt", sep = "\t")
peptides = pd.read_csv("peptides.txt", sep = "\t")
protein = pd.read_csv("proteinGroups.csv", sep = "\t")

# stats

len(peptides_tryptic)
len(peptides)
len(protein_tryptic)
len(protein)

peptides_tryptic.columns
peptides.columns
protein_tryptic.columns
protein.columns


protein





for i in protein.columns[0:50]:
    print(i)
    
    
    
    
df = pd.read_csv("peptide_tryptic_melted.csv", sep = "\t")
























