#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:39:01 2021

@author: ptruong
"""

def get_cell_line_state_replicate():
    cell_lines = ["A549", "MCF7", "RKO"]
    states = ["S", "D"]
    replicates = [1,2,3]
    return cell_lines, states, replicates

def get_base_cols_peptide():
    cols = ["Sequence", "Missed cleavages", "Proteins", "Leading razor protein",
       "Unique (Groups)", "Unique (Proteins)", "PEP", "Score"] 
    return cols

def get_base_cols_proteinGroups():
    base_cols = ["Protein IDs", "Majority protein IDs", "Peptide counts (all)", 
             "Peptide counts (razor+unique)", "Peptide counts (unique)", "Protein names",
             "Number of proteins", "Peptides", "Razor + unique peptides", "Unique peptides",
             "Sequence coverage [%]", "Unique + razor sequence coverage [%]",
             "Unique sequence coverage [%]", "Q-value", "Score"]
    return base_cols


def get_all_peptide_counts():
    peptides_counts_col = ["Peptides A549_D_Rep1", "Peptides A549_D_Rep2", "Peptides A549_D_Rep3", 
    "Peptides A549_S_Rep1", "Peptides A549_S_Rep2", "Peptides A549_S_Rep3",
    "Peptides MCF7_D_Rep1","Peptides MCF7_D_Rep2","Peptides MCF7_D_Rep3",
    "Peptides MCF7_S_Rep1","Peptides MCF7_S_Rep2","Peptides MCF7_S_Rep3",
    "Peptides RKO_D_Rep1","Peptides RKO_D_Rep2","Peptides RKO_D_Rep3",
    "Peptides RKO_S_Rep1","Peptides RKO_S_Rep2","Peptides RKO_S_Rep3"]
    return peptides_counts_col

def get_razor_and_unique_peptide_counts():
    razor_and_unique_peptide_counts_col = ["Razor + unique peptides A549_D_Rep1","Razor + unique peptides A549_D_Rep2","Razor + unique peptides A549_D_Rep3",
                                      "Razor + unique peptides A549_S_Rep1","Razor + unique peptides A549_S_Rep2","Razor + unique peptides A549_S_Rep3",
                                      "Razor + unique peptides MCF7_D_Rep1","Razor + unique peptides MCF7_D_Rep2","Razor + unique peptides MCF7_D_Rep3",
                                      "Razor + unique peptides MCF7_S_Rep1","Razor + unique peptides MCF7_S_Rep2","Razor + unique peptides MCF7_S_Rep3",
                                      "Razor + unique peptides RKO_D_Rep1","Razor + unique peptides RKO_D_Rep2","Razor + unique peptides RKO_D_Rep3",
                                      "Razor + unique peptides RKO_S_Rep1","Razor + unique peptides RKO_S_Rep2","Razor + unique peptides RKO_S_Rep3"]
    return razor_and_unique_peptide_counts_col

def get_unique_peptides():
    unique_peptides = ["Unique peptides A549_D_Rep1","Unique peptides A549_D_Rep2","Unique peptides A549_D_Rep3",
                       "Unique peptides A549_S_Rep1","Unique peptides A549_S_Rep2","Unique peptides A549_S_Rep3",
                       "Unique peptides MCF7_D_Rep1","Unique peptides MCF7_D_Rep2","Unique peptides MCF7_D_Rep3",
                       "Unique peptides MCF7_S_Rep1","Unique peptides MCF7_S_Rep2","Unique peptides MCF7_S_Rep3",
                       "Unique peptides RKO_D_Rep1","Unique peptides RKO_D_Rep2","Unique peptides RKO_D_Rep3",
                       "Unique peptides RKO_S_Rep1","Unique peptides RKO_S_Rep2","Unique peptides RKO_S_Rep3"]
    return unique_peptides

def get_sequence_coverage():
    sequence_coverage_cols = ["Sequence coverage A549_D_Rep1 [%]","Sequence coverage A549_D_Rep2 [%]","Sequence coverage A549_D_Rep3 [%]",
    "Sequence coverage A549_S_Rep1 [%]","Sequence coverage A549_S_Rep2 [%]","Sequence coverage A549_S_Rep3 [%]",
    "Sequence coverage MCF7_D_Rep1 [%]","Sequence coverage MCF7_D_Rep2 [%]","Sequence coverage MCF7_D_Rep3 [%]",
    "Sequence coverage MCF7_S_Rep1 [%]","Sequence coverage MCF7_S_Rep2 [%]","Sequence coverage MCF7_S_Rep3 [%]",
    "Sequence coverage RKO_D_Rep1 [%]","Sequence coverage RKO_D_Rep2 [%]","Sequence coverage RKO_D_Rep3 [%]",
    "Sequence coverage RKO_S_Rep1 [%]","Sequence coverage RKO_S_Rep2 [%]","Sequence coverage RKO_S_Rep3 [%]"]
    return sequence_coverage_cols

def get_all_reporter_intensity_correct():
    cols = ['Reporter intensity corrected 0 A549_D_Rep1',
 'Reporter intensity corrected 1 A549_D_Rep1',
 'Reporter intensity corrected 2 A549_D_Rep1',
 'Reporter intensity corrected 3 A549_D_Rep1',
 'Reporter intensity corrected 4 A549_D_Rep1',
 'Reporter intensity corrected 5 A549_D_Rep1',
 'Reporter intensity corrected 6 A549_D_Rep1',
 'Reporter intensity corrected 7 A549_D_Rep1',
 'Reporter intensity corrected 8 A549_D_Rep1',
 'Reporter intensity corrected 9 A549_D_Rep1',
 'Reporter intensity corrected 0 A549_D_Rep2',
 'Reporter intensity corrected 1 A549_D_Rep2',
 'Reporter intensity corrected 2 A549_D_Rep2',
 'Reporter intensity corrected 3 A549_D_Rep2',
 'Reporter intensity corrected 4 A549_D_Rep2',
 'Reporter intensity corrected 5 A549_D_Rep2',
 'Reporter intensity corrected 6 A549_D_Rep2',
 'Reporter intensity corrected 7 A549_D_Rep2',
 'Reporter intensity corrected 8 A549_D_Rep2',
 'Reporter intensity corrected 9 A549_D_Rep2',
 'Reporter intensity corrected 0 A549_D_Rep3',
 'Reporter intensity corrected 1 A549_D_Rep3',
 'Reporter intensity corrected 2 A549_D_Rep3',
 'Reporter intensity corrected 3 A549_D_Rep3',
 'Reporter intensity corrected 4 A549_D_Rep3',
 'Reporter intensity corrected 5 A549_D_Rep3',
 'Reporter intensity corrected 6 A549_D_Rep3',
 'Reporter intensity corrected 7 A549_D_Rep3',
 'Reporter intensity corrected 8 A549_D_Rep3',
 'Reporter intensity corrected 9 A549_D_Rep3',
 'Reporter intensity corrected 0 A549_S_Rep1',
 'Reporter intensity corrected 1 A549_S_Rep1',
 'Reporter intensity corrected 2 A549_S_Rep1',
 'Reporter intensity corrected 3 A549_S_Rep1',
 'Reporter intensity corrected 4 A549_S_Rep1',
 'Reporter intensity corrected 5 A549_S_Rep1',
 'Reporter intensity corrected 6 A549_S_Rep1',
 'Reporter intensity corrected 7 A549_S_Rep1',
 'Reporter intensity corrected 8 A549_S_Rep1',
 'Reporter intensity corrected 9 A549_S_Rep1',
 'Reporter intensity corrected 0 A549_S_Rep2',
 'Reporter intensity corrected 1 A549_S_Rep2',
 'Reporter intensity corrected 2 A549_S_Rep2',
 'Reporter intensity corrected 3 A549_S_Rep2',
 'Reporter intensity corrected 4 A549_S_Rep2',
 'Reporter intensity corrected 5 A549_S_Rep2',
 'Reporter intensity corrected 6 A549_S_Rep2',
 'Reporter intensity corrected 7 A549_S_Rep2',
 'Reporter intensity corrected 8 A549_S_Rep2',
 'Reporter intensity corrected 9 A549_S_Rep2',
 'Reporter intensity corrected 0 A549_S_Rep3',
 'Reporter intensity corrected 1 A549_S_Rep3',
 'Reporter intensity corrected 2 A549_S_Rep3',
 'Reporter intensity corrected 3 A549_S_Rep3',
 'Reporter intensity corrected 4 A549_S_Rep3',
 'Reporter intensity corrected 5 A549_S_Rep3',
 'Reporter intensity corrected 6 A549_S_Rep3',
 'Reporter intensity corrected 7 A549_S_Rep3',
 'Reporter intensity corrected 8 A549_S_Rep3',
 'Reporter intensity corrected 9 A549_S_Rep3',
 'Reporter intensity corrected 0 MCF7_D_Rep1',
 'Reporter intensity corrected 1 MCF7_D_Rep1',
 'Reporter intensity corrected 2 MCF7_D_Rep1',
 'Reporter intensity corrected 3 MCF7_D_Rep1',
 'Reporter intensity corrected 4 MCF7_D_Rep1',
 'Reporter intensity corrected 5 MCF7_D_Rep1',
 'Reporter intensity corrected 6 MCF7_D_Rep1',
 'Reporter intensity corrected 7 MCF7_D_Rep1',
 'Reporter intensity corrected 8 MCF7_D_Rep1',
 'Reporter intensity corrected 9 MCF7_D_Rep1',
 'Reporter intensity corrected 0 MCF7_D_Rep2',
 'Reporter intensity corrected 1 MCF7_D_Rep2',
 'Reporter intensity corrected 2 MCF7_D_Rep2',
 'Reporter intensity corrected 3 MCF7_D_Rep2',
 'Reporter intensity corrected 4 MCF7_D_Rep2',
 'Reporter intensity corrected 5 MCF7_D_Rep2',
 'Reporter intensity corrected 6 MCF7_D_Rep2',
 'Reporter intensity corrected 7 MCF7_D_Rep2',
 'Reporter intensity corrected 8 MCF7_D_Rep2',
 'Reporter intensity corrected 9 MCF7_D_Rep2',
 'Reporter intensity corrected 0 MCF7_D_Rep3',
 'Reporter intensity corrected 1 MCF7_D_Rep3',
 'Reporter intensity corrected 2 MCF7_D_Rep3',
 'Reporter intensity corrected 3 MCF7_D_Rep3',
 'Reporter intensity corrected 4 MCF7_D_Rep3',
 'Reporter intensity corrected 5 MCF7_D_Rep3',
 'Reporter intensity corrected 6 MCF7_D_Rep3',
 'Reporter intensity corrected 7 MCF7_D_Rep3',
 'Reporter intensity corrected 8 MCF7_D_Rep3',
 'Reporter intensity corrected 9 MCF7_D_Rep3',
 'Reporter intensity corrected 0 MCF7_S_Rep1',
 'Reporter intensity corrected 1 MCF7_S_Rep1',
 'Reporter intensity corrected 2 MCF7_S_Rep1',
 'Reporter intensity corrected 3 MCF7_S_Rep1',
 'Reporter intensity corrected 4 MCF7_S_Rep1',
 'Reporter intensity corrected 5 MCF7_S_Rep1',
 'Reporter intensity corrected 6 MCF7_S_Rep1',
 'Reporter intensity corrected 7 MCF7_S_Rep1',
 'Reporter intensity corrected 8 MCF7_S_Rep1',
 'Reporter intensity corrected 9 MCF7_S_Rep1',
 'Reporter intensity corrected 0 MCF7_S_Rep2',
 'Reporter intensity corrected 1 MCF7_S_Rep2',
 'Reporter intensity corrected 2 MCF7_S_Rep2',
 'Reporter intensity corrected 3 MCF7_S_Rep2',
 'Reporter intensity corrected 4 MCF7_S_Rep2',
 'Reporter intensity corrected 5 MCF7_S_Rep2',
 'Reporter intensity corrected 6 MCF7_S_Rep2',
 'Reporter intensity corrected 7 MCF7_S_Rep2',
 'Reporter intensity corrected 8 MCF7_S_Rep2',
 'Reporter intensity corrected 9 MCF7_S_Rep2',
 'Reporter intensity corrected 0 MCF7_S_Rep3',
 'Reporter intensity corrected 1 MCF7_S_Rep3',
 'Reporter intensity corrected 2 MCF7_S_Rep3',
 'Reporter intensity corrected 3 MCF7_S_Rep3',
 'Reporter intensity corrected 4 MCF7_S_Rep3',
 'Reporter intensity corrected 5 MCF7_S_Rep3',
 'Reporter intensity corrected 6 MCF7_S_Rep3',
 'Reporter intensity corrected 7 MCF7_S_Rep3',
 'Reporter intensity corrected 8 MCF7_S_Rep3',
 'Reporter intensity corrected 9 MCF7_S_Rep3',
 'Reporter intensity corrected 0 RKO_D_Rep1',
 'Reporter intensity corrected 1 RKO_D_Rep1',
 'Reporter intensity corrected 2 RKO_D_Rep1',
 'Reporter intensity corrected 3 RKO_D_Rep1',
 'Reporter intensity corrected 4 RKO_D_Rep1',
 'Reporter intensity corrected 5 RKO_D_Rep1',
 'Reporter intensity corrected 6 RKO_D_Rep1',
 'Reporter intensity corrected 7 RKO_D_Rep1',
 'Reporter intensity corrected 8 RKO_D_Rep1',
 'Reporter intensity corrected 9 RKO_D_Rep1',
 'Reporter intensity corrected 0 RKO_D_Rep2',
 'Reporter intensity corrected 1 RKO_D_Rep2',
 'Reporter intensity corrected 2 RKO_D_Rep2',
 'Reporter intensity corrected 3 RKO_D_Rep2',
 'Reporter intensity corrected 4 RKO_D_Rep2',
 'Reporter intensity corrected 5 RKO_D_Rep2',
 'Reporter intensity corrected 6 RKO_D_Rep2',
 'Reporter intensity corrected 7 RKO_D_Rep2',
 'Reporter intensity corrected 8 RKO_D_Rep2',
 'Reporter intensity corrected 9 RKO_D_Rep2',
 'Reporter intensity corrected 0 RKO_D_Rep3',
 'Reporter intensity corrected 1 RKO_D_Rep3',
 'Reporter intensity corrected 2 RKO_D_Rep3',
 'Reporter intensity corrected 3 RKO_D_Rep3',
 'Reporter intensity corrected 4 RKO_D_Rep3',
 'Reporter intensity corrected 5 RKO_D_Rep3',
 'Reporter intensity corrected 6 RKO_D_Rep3',
 'Reporter intensity corrected 7 RKO_D_Rep3',
 'Reporter intensity corrected 8 RKO_D_Rep3',
 'Reporter intensity corrected 9 RKO_D_Rep3',
 'Reporter intensity corrected 0 RKO_S_Rep1',
 'Reporter intensity corrected 1 RKO_S_Rep1',
 'Reporter intensity corrected 2 RKO_S_Rep1',
 'Reporter intensity corrected 3 RKO_S_Rep1',
 'Reporter intensity corrected 4 RKO_S_Rep1',
 'Reporter intensity corrected 5 RKO_S_Rep1',
 'Reporter intensity corrected 6 RKO_S_Rep1',
 'Reporter intensity corrected 7 RKO_S_Rep1',
 'Reporter intensity corrected 8 RKO_S_Rep1',
 'Reporter intensity corrected 9 RKO_S_Rep1',
 'Reporter intensity corrected 0 RKO_S_Rep2',
 'Reporter intensity corrected 1 RKO_S_Rep2',
 'Reporter intensity corrected 2 RKO_S_Rep2',
 'Reporter intensity corrected 3 RKO_S_Rep2',
 'Reporter intensity corrected 4 RKO_S_Rep2',
 'Reporter intensity corrected 5 RKO_S_Rep2',
 'Reporter intensity corrected 6 RKO_S_Rep2',
 'Reporter intensity corrected 7 RKO_S_Rep2',
 'Reporter intensity corrected 8 RKO_S_Rep2',
 'Reporter intensity corrected 9 RKO_S_Rep2',
 'Reporter intensity corrected 0 RKO_S_Rep3',
 'Reporter intensity corrected 1 RKO_S_Rep3',
 'Reporter intensity corrected 2 RKO_S_Rep3',
 'Reporter intensity corrected 3 RKO_S_Rep3',
 'Reporter intensity corrected 4 RKO_S_Rep3',
 'Reporter intensity corrected 5 RKO_S_Rep3',
 'Reporter intensity corrected 6 RKO_S_Rep3',
 'Reporter intensity corrected 7 RKO_S_Rep3',
 'Reporter intensity corrected 8 RKO_S_Rep3',
 'Reporter intensity corrected 9 RKO_S_Rep3']
    return cols

def get_reporter_intensity_without_control():
    cols = get_all_reporter_intensity_correct()
    cols_without_control = []
    for i in cols:
        if i.split(" ")[3] != "0":
            cols_without_control.append(i)
    return cols_without_control
