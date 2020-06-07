#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:48:58 2020

@author: ptruong
"""

import pandas as pd
import numpy as np
from helper import *



df = pd.read_csv("proteinGroups_tryptic_head2000.csv", sep = "\t")

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
# df = df[df["PEP"] < 0.01]


# Razor peptide?
# https://med.uottawa.ca/core-facilities/facilities/proteomic/resources/interpret-result


# Q-value treshold - protein
df = df[df["Q-value"] < 0.01]

# Peptides treshold - protein
df = df[df["Peptides"] > 1]

##########################################
# map peptides to report intensities NaN #
##########################################
treshold = 1

if df["Unique peptides A549_D_Rep1"] > treshold:
    df["Reporter intensity corrected 0 A549_D_Rep1"]

df["Unique peptides A549_D_Rep2"] > treshold
df["Unique peptides A549_D_Rep3"] > treshold
df["Unique peptides A549_S_Rep1"] > treshold
df["Unique peptides A549_S_Rep2"] > treshold
df["Unique peptides A549_S_Rep3"] > treshold
df["Unique peptides MCF7_D_Rep1"] > treshold
df["Unique peptides MCF7_D_Rep2"] > treshold
df["Unique peptides MCF7_D_Rep3"] > treshold
df["Unique peptides MCF7_S_Rep1"] > treshold
df["Unique peptides MCF7_S_Rep2"] > treshold
df["Unique peptides MCF7_S_Rep3"] > treshold
df["Unique peptides RKO_D_Rep1"] > treshold
df["Unique peptides RKO_D_Rep2"] > treshold
df["Unique peptides RKO_D_Rep3"] > treshold
df["Unique peptides RKO_S_Rep1"] > treshold
df["Unique peptides RKO_S_Rep2"] > treshold
df["Unique peptides RKO_S_Rep3"] > treshold

df["Reporter intensity corrected 0 A549_D_Rep1"] * (df["Unique peptides A549_D_Rep1"] > treshold)
(df["Unique peptides A549_D_Rep1"] > treshold)
df["Reporter intensity corrected 0 A549_D_Rep1"]


def treshold_by_peptides(intensities, peptide_counts, treshold = 1):
    #intensities = df["Reporter intensity corrected 0 A549_D_Rep1"]
    #booleans = (df["Unique peptides A549_D_Rep1"] > treshold)
    
    intensities = intensities
    peptide_counts = peptide_counts
    treshold = treshold
    
    booleans = (peptide_counts > treshold)
    res = intensities * booleans.apply(lambda x: np.nan if x == False else True)
    return res

treshold_by_peptides(df["Reporter intensity corrected 0 A549_D_Rep1"], df["Unique peptides A549_D_Rep1"])

peptide_count_1 = df["Peptides A549_D_Rep1"]
peptide_count_2 = df["Peptides A549_D_Rep2"]

intensities_1 = df["Reporter intensity corrected 0 A549_D_Rep1"]
intensities_2 = df["Reporter intensity corrected 0 A549_D_Rep2"]

t_1 = treshold_by_peptides(intensities_1, peptide_count_1)
t_2 = treshold_by_peptides(intensities_2, peptide_count_2)


def gen_cell_lines_states_replicates():
    cell_lines = ["A549", "MCF7", "RKO"]
    states = ["D", "S"]
    replicates = ["Rep1", "Rep2", "Rep3"]
    var_list = []
    for i in range(len(cell_lines)):
        for j in range(len(states)):
            for k in range(len(replicates)):
                var_str = ""
                var_str += cell_lines[i] + "_"
                var_str += states[j] + "_"
                var_str += replicates[k]
                var_list.append(var_str)
    return var_list

def add_prefix_with_treatments(prefix = "Reporter intensity corrected", treatments = 10):
    var_list = gen_cell_lines_states_replicates()
    treatments = [i for i in range(treatments)]
    prefix = prefix
    res_list = []
    for i in var_list:
        for j in treatments:
            unit_str = prefix + " "
            unit_str += str(j) + " "
            unit_str += i
            res_list.append(unit_str)
    return res_list

def add_prefix(prefix = "Peptides"):
    var_list = gen_cell_lines_states_replicates()
    prefix = prefix
    res_list = []
    for i in var_list:
        unit_str = prefix + " "
        unit_str += i
        res_list.append(unit_str)
    return res_list

reporter_intensity_corrected = add_prefix_with_treatments()
unique_peptides = add_prefix()

data = pd.DataFrame()
#data[df[i].name] = df[i].values
#data[df["Reporter intensity corrected 9 RKO_S_Rep2"].name] = df["Reporter intensity corrected 9 RKO_S_Rep2"].values
"""
for i in reporter_intensity_corrected[:]:
    for j in unique_peptides[:]:
        print(i + "####" + j)
        data[df[i].name] = treshold_by_peptides(df[i], df[j])
"""

for i in unique_peptides:
    for j in reporter_intensity_corrected:
        if i.split()[-1] == j.split()[-1]:
            print(i + "\t - \t" + j)
            data[df[j].name] = treshold_by_peptides(df[j], df[i])
            
#Count Na --> 

print(df.head())

peptide_count_1 = df["Peptides A549_D_Rep1"]
peptide_count_2 = df["Peptides A549_D_Rep2"]

intensities_1 = df["Reporter intensity corrected 0 A549_D_Rep1"]
intensities_2 = df["Reporter intensity corrected 0 A549_D_Rep2"]

t_1 = treshold_by_peptides(intensities_1, peptide_count_1)
t_2 = treshold_by_peptides(intensities_2, peptide_count_2)

def count_nan(df):
    count_nan = len(df) - df.count()
    return count_nan

treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])
treshold_by_peptides(df["Reporter intensity corrected "])




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

# Check to check tresholding
df.iloc[:,13:20]



df["Unique peptides A549_D_Rep1"]
df["Razor + unique peptides A549_D_Rep1"]
df["Reporter intensity corrected 0 A549_D_Rep1"]











