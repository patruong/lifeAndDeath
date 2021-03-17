#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:16:03 2021

@author: ptruong
"""

import pandas as pd
import numpy as np

from get_columns import get_cell_line_state_replicate, get_base_cols_proteinGroups, get_all_peptide_counts, get_razor_and_unique_peptide_counts, get_unique_peptides, get_sequence_coverage, get_all_reporter_intensity_correct

df = pd.read_csv("proteinGroups tryptic.csv", sep = "\t")

cell_lines, states, replicates = get_cell_line_state_replicate()

base_cols = get_base_cols_proteinGroups()
peptide_count_cols = get_all_peptide_counts()
#razor_and_unique_peptides_cols = get_razor_and_unique_peptide_counts()
#unique_peptides_cols = get_unique_peptides()
sequence_coverage_cols = get_sequence_coverage()
reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()

def select_rep_state_intensities(rep, state):
    intensity_list = []
    for i in reporter_intensity_corrected_cols:
        if i.split(" ")[4].split("_")[2] == ("Rep" + str(rep)):
            if i.split(" ")[4].split("_")[1] == state:
                intensity_list.append(i)
    return intensity_list

def select_rep_state_cell_line_intensities(rep, state, cell_line):
    intensity_list = []
    for i in reporter_intensity_corrected_cols:
        if i.split(" ")[4].split("_")[2] == ("Rep" + str(rep)):
            if i.split(" ")[4].split("_")[1] == state:
                if i.split(" ")[4].split("_")[0] == cell_line:
                    intensity_list.append(i)
    return intensity_list

df_base = df[get_base_cols_proteinGroups()]
df_peptide_count = df[peptide_count_cols]
df_sequence_coverage = df[sequence_coverage_cols]


def apply_peptide_count_treshold(df_subset, treshold = 1):
    for i in df_subset.columns:
        col = "Peptides " + i.split(" ")[-1]
        peptide_count = df_peptide_count[col]
        peptide_count_treshold_boolean_array = (peptide_count > treshold)
        df_subset[i] = df_subset[i] * peptide_count_treshold_boolean_array
    return df_subset

def apply_sequence_coverage_treshold(df_subset, treshold = 0):
    for i in df_subset.columns:
        col = "Sequence coverage " + i.split(" ")[-1] + " [%]"
        sequence_coverage = df_sequence_coverage[col]
        sequence_coverage_treshold_boolean_array = (sequence_coverage > treshold)
        df_subset[i] = df_subset[i] * sequence_coverage_treshold_boolean_array
    return df_subset

def apply_treshold(df_subset, peptide_count_treshold = 1, sequence_coverage_percentage_treshold = 0):
    df_subset = apply_peptide_count_treshold(df_subset, peptide_count_treshold)
    df_subset = apply_sequence_coverage_treshold(df_subset, sequence_coverage_percentage_treshold)
    return df_subset


peptide_count_treshold = 1
sequence_coverage_percentage_treshold = 0
df_t = apply_treshold(df[reporter_intensity_corrected_cols], peptide_count_treshold, sequence_coverage_percentage_treshold)
df_t = df_t.replace({0:np.nan})
df_t = np.log2(df_t)


df_subset = df_t[select_rep_state_cell_line_intensities(1, "D", "RKO")]

control = pd.DataFrame(df_subset[df_subset.columns[0]])

df_fc = df_subset - control









