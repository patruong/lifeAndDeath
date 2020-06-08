#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 15:11:58 2020

@author: ptruong
"""

import pandas as pd
import numpy as np
import os 

files = os.listdir()
print(os.listdir())
files[1] # peptides tryptic.txt
files[2] 
files[3] # proteinGroups tryptic.txt
files[4] # proteinGroups.txt
files[5] # peptides.txt


# Using readlines() 
file1 = open(files[1], 'r') 
lines = file1.readline() 
print(lines)
print(lines.split("\t"))
print(len(lines.split("\t")))
df.index



df = pd.read_csv(files[1], sep = "\t")



df["Gene names"]

df["PEP"] # Need to PEP treshold

# Need to normalize

# PEP treshold
df = df[df["PEP"] < 0.01]

# Normalize
df['Identification type A549_D_Rep1'

start = 680
stop = start+20
for i in df.columns[start:stop]:
    print(i)


df['Reporter intensity 9 A549_D_Rep3']
df['Reporter intensity corrected 9 A549_D_Rep3']
df['Reporter intensity count 9 A549_D_Rep3']
df['Reverse'].unique()
df['Potential contaminant'].unique()
df['id'].unique()
df['Protein group IDs'].unique()

# Lets try Corrected
df['Reporter intensity corrected 9 A549_D_Rep3']


df.drop(df.iloc[:, 1:98], inplace = True, axis = 1) 


df["Experiment A549_D_Rep1"]
df['Reporter intensity corrected 0']


##################
# Generates cols #
##################

def pick_cols(init_string = "Reporter intensity corrected ", drug_codes = [i for i in range(10)], cell_lines = ["A549", "MCF7", "RKO"], states = ["D", "S"], replicates = ["Rep{}".format(i+1) for i in range(3)]):
    # 0 - 9 , TMT and the drugs
    drug_codes = drug_codes
    n_drug_codes = len(drug_codes)
    
    # cell_lines
    cell_lines = cell_lines
    n_cell_lines = len(cell_lines)
    
    # states
    states = states
    n_states = len(states)
    
    # replicates
    replicates = replicates
    n_replicates = len(replicates)
    
    # Columns
    #init_strings = ["Reporter intensity ", "Reporter intensity corrected ", "Reporter intensity count "]
    #init_string = init_strings[init_string_number]
    init_string = init_string
    cols = []
    for i in range(n_drug_codes):
        for j in range(n_cell_lines):
            for k in range(n_states):
                for l in range(n_replicates):
                    col_str = init_string
                    col_str += str(drug_codes[i]) + " "
                    col_str += str(cell_lines[j]) + "_"
                    col_str += states[k] + "_"
                    col_str += replicates[l]
                    cols.append(col_str)
    return cols
       
        
cols = generate_cols(1)
df_ = df[cols]

df_.to_csv("testfile.csv", sep = ";")
       
        
        


