#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:02:28 2021

@author: ptruong
"""



import pandas as pd 
import numpy as np
from get_columns import get_cell_line_state_replicate, get_base_cols_peptide, get_all_reporter_intensity_correct, get_reporter_intensity_without_control()

df = pd.read_csv("peptides tryptic.csv", sep = "\t")

cell_lines, states, replicates = get_cell_line_state_replicate()
base_cols = get_base_cols_peptide()
reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()

df_base = df[get_base_cols_peptide()]


def select_rep_state_cell_line_intensities(rep, state, cell_line):
    intensity_list = []
    for i in reporter_intensity_corrected_cols:
        if i.split(" ")[4].split("_")[2] == ("Rep" + str(rep)):
            if i.split(" ")[4].split("_")[1] == state:
                if i.split(" ")[4].split("_")[0] == cell_line:
                    intensity_list.append(i)
    return intensity_list


def subtract_df_with_col(df, col):
    """
    test to see results:
    a = pd.DataFrame(np.array([[1,1,1],[2,4,8],[3,6,9]]).T, columns = ["a", "b", "c"])
    (a.T - a.a).T
    subtract_df_with_col(a, "a")
    """
    return (df.T - df[col]).T

def compute_control_to_treated_fc(df_subset):
    control_col = df_subset.columns[0]
    fc = subtract_df_with_col(df_subset, control_col)
    return fc


def get_df_with_control_vs_treated_fc(df_t):
    """
    df_t is tresholded and logged df with all reporter_intensity_corrected_cols
    """
    df_fc = pd.DataFrame()
    for cell_line in cell_lines:
        for state in states:
            for replicate in replicates:
                df_subset = df_t[select_rep_state_cell_line_intensities(replicate, state, cell_line)]            
                fc = compute_control_to_treated_fc(df_subset)
                if df_fc.empty:
                    df_fc = df_fc.append(fc)
                else:
                    df_fc = df_fc.join(fc)
    return df_fc

df_int = df[reporter_intensity_corrected_cols]
df_int = df_int.replace({0:np.nan})
df_int = np.log2(df_int)
df_fc = get_df_with_control_vs_treated_fc(df_int) #control - treated, therefore control channels are 0
df_res = df_base.join(df_fc)


#peptides without treshold:225715
# PEP treshold
df_res = df_res[df_res["PEP"] < 0.05] # We needed to do this step to treshold for q-values. perhaps could have done earlier?
# peptides after q-value treshold: 216137
# Misssed cleavages treshold
df_res = df_res[df_res["Missed cleavages"] < 2]
# peptides after missed cleavages treshold: 210031

# This is the data matrix we work with.
df_int = df_res[get_reporter_intensity_without_control()]









