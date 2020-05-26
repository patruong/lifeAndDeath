#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:48:58 2020

@author: ptruong
"""

import pandas as pd
import numpy as np
from helper import *



df = pd.read_csv("proteinGroups tryptic.txt", sep = "\t")

# peptide file
df["MS/MS Count"]
df["Best MS/MS"]

# protein file
df["Peptides"]
df["Peptide counts (all)"].unique()
df["Fraction 1"].unique()
df["Sequence coverage [%]"].min()
df["Razor + unique peptides A549_D_Rep1"].unique()
df["Unique peptides A549_D_Rep1"].unique()


start = 0
stop = start + 200
for i in df.columns[start:stop]:
    try:
        print(i)
    except:
        print("END")
        break
    
# PEP treshold - peptide
#df = df[df["PEP"] < 0.01]


# Razor peptide?
# https://med.uottawa.ca/core-facilities/facilities/proteomic/resources/interpret-result


# Q-value treshold - protein
df = df[df["Q-value"] < 0.01]

# Peptides treshold - protein
df = df[df["Peptides"] > 1]


# pick relevant cols
cols = pick_cols(init_string = "Reporter intensity corrected ", 
                 drug_codes = [i for i in range(0, 10)], 
                 cell_lines = ["A549", "MCF7", "RKO"], 
                 states = ["D", "S"], 
                 replicates = ["Rep{}".format(i+1) for i in range(3)])
df = df[cols]

# Normalize ( mean normalization because outliers might screw up 0-1 normalization)
normalized_df=(df-df.mean())/df.std()

# Take subset



df["Unique peptides A549_D_Rep1"]
df["Razor + unique peptides A549_D_Rep1"]
df["Reporter intensity corrected 0 A549_D_Rep1"]












