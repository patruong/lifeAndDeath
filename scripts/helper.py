#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:06:59 2020

@author: ptruong
"""

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
       