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


def compute_log2FC(df_int):
    cell_lines = [ "A549", "MCF7", "RKO"]
    states = ["S", "D"]
    replicates = ["Rep1", "Rep2", "Rep3"]
    treatments = [i for i in range(10)]

    df_int = df_int.replace(0, np.nan)
    df_log2 = np.log2(df_int)

    fc_array = []
    for cell_line in cell_lines:
        for state in states:
            for rep in replicates:
                df_sample = df_log2.iloc[:,df_log2.columns.get_level_values("batch") == f"{cell_line}_{state}_{rep}"]
                case = df_sample.iloc[:, df_sample.columns.get_level_values("treatment") != "0"]
                control = df_sample.iloc[:, df_sample.columns.get_level_values("treatment") == "0"]
                fc = pd.DataFrame((control.values - case.values), index = case.index, columns = case.columns)
                fc_array.append(fc)


    df_fc = pd.concat(fc_array, axis = 1)
    return df_fc

def split_surviving_dead(df_fc, cell_lines = ["A549", "MCF7", "RKO"]):
    # df_fc = df_log2FC
    S = df_fc[cell_lines].iloc[:, (df_fc[cell_lines].columns.get_level_values("state") == "S")]
    D = df_fc[cell_lines].iloc[:, df_fc[cell_lines].columns.get_level_values("state") == "D"]
    S = pd.DataFrame(S.values, index = S.index, columns = S.columns.get_level_values("experiment"))
    D = pd.DataFrame(D.values, index = D.index, columns = D.columns.get_level_values("experiment"))
    col_mapper = lambda x: x.split("_")[0] + "_" + x.split("_")[1] + "_" + x.split("_")[-1] 
    S = S.rename(columns = col_mapper)
    D = D.rename(columns = col_mapper)
    SD = pd.concat([S.stack() ,D.stack()], axis = 1).rename(columns={0:"S", 1:"D"})
    return S, D, SD

