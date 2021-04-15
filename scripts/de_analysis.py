#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 21:49:28 2021

@author: ptruong
"""

import numpy as np
import pandas as pd

from midx import col_to_mIdx, intensities_to_midx_df, diffacto_col_to_mIdx, diffacto_to_midx_df
from top3 import top3, protSum_col_to_mIdx, protSum_intensities_to_midx_df, aggregate_protein_quantity, get_p_matrix

def get_variables():
    cell_lines = ['A549', 'MCF7', 'RKO']
    states = ['D', 'S']
    replicates = ['Rep1', 'Rep2', 'Rep3']
    return cell_lines, states, replicates


def get_log2FC_regulation(df_quant):
    cell_lines, states, replicates = get_variables()
    df_regulation = pd.DataFrame()
    for cell_line in cell_lines:
        for state in states:
            for treatment_i in range(1,10):
                df_reg_sample = pd.DataFrame(np.log2(df_quant[cell_line][state][str(treatment_i)]) - np.log2(df_quant[cell_line][state]["0"]), 
                                             columns = pd.MultiIndex(levels = [[cell_line], [state], [str(treatment_i)]], 
                                                                     codes = [[0], [0], [0]], 
                                                                     name = ["cell_line", "state", "treatment"]))
                df_regulation = pd.concat([df_regulation, df_reg_sample], axis = 1)
    return df_regulation


def get_log2FC(df_protQuant):
    df_quant = aggregate_protein_quantity(df_protQuant)
    df_quant = diffacto_to_midx_df(df_quant, protein_index=False)
    df_pVals, df_tStats = get_p_matrix(df_protQuant)
    df_log2FC = get_log2FC_regulation(df_quant)
    df_pVals, df_tStats = diffacto_to_midx_df(df_pVals, protein_index=False), diffacto_to_midx_df(df_tStats, protein_index=False)
    return df_log2FC, df_pVals, df_tStats
