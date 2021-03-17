#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:02:28 2021

@author: ptruong
"""


from functools import reduce

import pandas as pd 
from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

df = pd.read_csv("peptides tryptic.csv", sep = "\t")

rep_list = list(pd.read_csv("reporter_intensity_list.csv", sep = "\t").columns)

def get_base_cols():
    cols = ["Sequence", "Missed cleavages", "Proteins", "Leading razor protein",
       "Unique (Groups)", "Unique (Proteins)", "PEP", "Score"] 
    return cols

def get_all_intensities():
    rep_list_corrected = []
    for i in rep_list:
        if i.split(" ")[2] == "corrected":
            rep_list_corrected.append(i)
    return rep_list_corrected

def select_rep_state_intensities(rep, state):
    rep = rep
    state = state
    
    rep_list_corrected = []
    for i in rep_list:
        if i.split(" ")[2] == "corrected":
            if i.split(" ")[4].split("_")[2] == ("Rep" + str(rep)):
                if i.split(" ")[4].split("_")[1] == state:
                    rep_list_corrected.append(i)
    return rep_list_corrected

df_base = df[get_base_cols()]
df_base = df

s1 = df[select_rep_state_intensities(1, "S")]
s2 = df[select_rep_state_intensities(2, "S")]
s3 = df[select_rep_state_intensities(3, "S")]
d1 = df[select_rep_state_intensities(1, "D")]
d2 = df[select_rep_state_intensities(2, "D")]
d3 = df[select_rep_state_intensities(3, "D")]


