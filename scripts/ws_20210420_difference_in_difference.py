#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 17:47:40 2021

@author: ptruong
"""



"""
Abstract:
    
    This script is to analyze the log2FC for unprocessed data.
    The theory is that since we are looking within each channel we should
    not need to normalize the data. 
    
    We are looking for proteins that behave different for different treatment.
    In other wrods, proteins that are up-regulated in one treatment, but
    down for another. It should behave similarily across technical replicates.
    We are in other words looking for difference in difference.
    
    
    We are looking into using MANOVA for this.
"""


import os
import time

import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from load import load_peptide_data
from top3 import top3
from get_variables import get_cell_lines_states_replicates

os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")
df = load_peptide_data("peptides tryptic.csv", max_PEP = 0.1, max_missed_cleavages = 3)

df = df.iloc[df.index.get_level_values("Unique (Groups)") == 'yes', :] # Take only peptides with unique protein groups

cell_lines,states, replicates = get_cell_lines_states_replicates()
treatments = [str(i) for i in range(10)]

# new top3 - take only the top scoring peptide like triqler

start = time.time()
df_best_scoring_proteins = pd.DataFrame()
for protein in df.index.get_level_values("Leading razor protein").unique():
    df_protein = df.iloc[df.index.get_level_values("Leading razor protein") == protein, :]
    df_best_score = df_protein.iloc[np.where(df_protein.index.get_level_values("Score") == df_protein.index.get_level_values("Score").max())[0][0], :]
    df_best_score = pd.DataFrame(df_best_score).T
    df_best_score.index.set_names(df_protein.index.names, inplace = True)
    df_best_scoring_proteins = pd.concat([df_best_scoring_proteins, df_best_score], axis = 0)
    print(time.time()-start)
end = time.time()

print(end-start)

cell_line = "A549"
state = "S"
control = "0"
treatment = "1"

df_best_scoring_proteins = df_best_scoring_proteins.replace(0,np.nan) # zero-score is assumed to be non-measurements.
df_log2_protein = np.log2(df_best_scoring_proteins)

df_log2FC = pd.DataFrame()
for cell_line in cell_lines:
    for state in states:
        control_sample = "_".join([cell_line, state, control])
        df_control = df_log2_protein.iloc[:, df_log2_protein.columns.get_level_values("sample") == control_sample]
        for treatment in treatments:
            treatment_sample = "_".join([cell_line, state, treatment]) 
        
            df_treatment = df_log2_protein.iloc[:, df_log2_protein.columns.get_level_values("sample") == treatment_sample]
            df_log2FC_part = pd.DataFrame(df_treatment.values - df_control.values, index = df_treatment.index, columns = df_treatment.columns)
            df_log2FC = pd.concat([df_log2FC, df_log2FC_part], axis = 1)

df_log2FC.iloc[:, df_log2FC.columns.get_level_values("sample") == "A549_S_1"]
df_log2FC.iloc[:, df_log2FC.columns.get_level_values("sample") == "A549_D_1"]

#df_means = df_best_scoring_proteins.T.groupby("sample").mean().T
# now we have only one protein per sample, we can perform log2FC....

import scipy.stats as stats


a =  np.array([[9.87, 9.03, 6.81],
              [7.18, 8.35, 7.00],
              [8.39, 7.58, 7.68],
              [7.45, 6.33, 9.35],
              [6.41, 7.10, 9.33],
              [8.00, 8.24, 8.44]])
a = pd.DataFrame(a)
b = np.array([[6.35, 7.30, 7.16],
              [6.65, 6.68, 7.63],
              [5.72, 7.73, 6.72],
              [7.01, 9.19, 7.41],
              [7.75, 7.87, 8.30],
              [6.90, 7.97, 6.97]])    
b = pd.DataFrame(b)
c = np.array([[3.31, 8.77, 1.01],
              [8.25, 3.24, 3.62],
              [6.32, 8.81, 5.19],
              [7.48, 8.83, 8.91],
              [8.59, 6.01, 6.07],
              [3.07, 9.72, 7.48]])
c = pd.DataFrame(c)
stats.f_oneway(a,b)

stats.f_oneway(a,b,c)
