#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:46:42 2021

@author: ptruong
"""

import numpy as np
import pandas as pd

def col_to_mIdx(df_int):
    treatment_mapper = lambda x : x.split(" ")[3]
    specie_mapper = lambda x : x.split(" ")[4].split("_")[0]
    state_mapper = lambda x : x.split(" ")[4].split("_")[1]
    replicate_mapper = lambda x : x.split(" ")[4].split("_")[2]
    batch_mapper = lambda x : x.split(" ")[-1]
    experiment_mapper = lambda x: x.split(" ")[3] + "_" + x.split(" ")[-1] 
    sample_mapper = lambda x: x.split(" ")[4].split("_")[0] + "_" + x.split(" ")[4].split("_")[1] + "_" +  x.split(" ")[3]
    
    def map_values_to_value_list(value_list, values):
        """
        e.g. test below lists to see results.
        value_list = ["a", "b", "c"]
        values = ["c","c","a"]
        """
        return [value_list.index(x) for x in values]
    
    specie_col = df_int.columns.map(specie_mapper)
    specie_list = list(specie_col.unique())
    specie_code = map_values_to_value_list(specie_list, specie_col)
    
    state_col = df_int.columns.map(state_mapper)
    state_list = list(state_col.unique())
    state_code = map_values_to_value_list(state_list, state_col)
    
    replicate_col = df_int.columns.map(replicate_mapper)
    replicate_list = list(replicate_col.unique())
    replicate_code = map_values_to_value_list(replicate_list, replicate_col)
    
    treatment_col = df_int.columns.map(treatment_mapper)
    treatment_list = list(treatment_col.unique())
    treatment_code = map_values_to_value_list(treatment_list, treatment_col)
    
    batch_col = df_int.columns.map(batch_mapper)
    batch_list = list(batch_col.unique())
    batch_code = map_values_to_value_list(batch_list, batch_col)
    
    sample_col = df_int.columns.map(sample_mapper)
    sample_list = list(sample_col.unique())
    sample_code = map_values_to_value_list(sample_list, sample_col)

    experiment_col = df_int.columns.map(experiment_mapper)
    experiment_list = list(experiment_col.unique())
    experiment_code = map_values_to_value_list(experiment_list, experiment_col)

    
    midx = pd.MultiIndex(levels=[specie_list, state_list, treatment_list, replicate_list, batch_list, sample_list, experiment_list], 
                         codes=[specie_code, state_code, treatment_code, replicate_code, batch_code, sample_code, experiment_code],
                         names=["cell_line", "state", "treatment", "replicate", "batch", "sample", "experiment"])
    return midx

def intensities_to_midx_df(df_int):
    """
    df_int = df[reporter_intensity_corrected_cols]

    example on selecting from levels on multiIndex:
        df_int.iloc[:, df_int.columns.get_level_values(level=1) == "S"]
        
    Example used for this build.
    
    Column to midx
    https://riptutorial.com/pandas/example/18696/how-to-change-standard-columns-to-multiindex
    
    Select columns in multiIdx
    https://stackoverflow.com/questions/25189575/pandas-dataframe-select-columns-in-multiindex
    
    # build hierarchical df
    
    df_ = pd.DataFrame(np.random.randn(2,3), columns=['a','b','c'])
    midx = pd.MultiIndex(levels=[['zero', 'one'], ['x','y']], codes=[[1,1,0,],[1,0,1,]])
    df_.columns = midx
    
    
    list1 = ['a','b','c','d'] # values
    list2 = ['d','a', "c", "d"] # value list
    print([list2.index(x) for x in list1 if x in list2])

    """
    midx = col_to_mIdx(df_int)
    df_int.columns = midx
    return df_int



def diffacto_col_to_mIdx(df):
    cell_line_mapper = lambda x: x.split("_")[0]
    state_mapper = lambda x: x.split("_")[1]
    treatment_mapper = lambda x: x.split("_")[2]
    
    def map_values_to_value_list(value_list, values):
        """
        e.g. test below lists to see results.
        value_list = ["a", "b", "c"]
        values = ["c","c","a"]
        """
        return [value_list.index(x) for x in values]
    
    cell_line_col = df.columns.map(cell_line_mapper)
    cell_line_list = list(cell_line_col.unique())
    cell_line_code = map_values_to_value_list(cell_line_list, cell_line_col)
    
    state_col = df.columns.map(state_mapper)
    state_list = list(state_col.unique())
    state_code = map_values_to_value_list(state_list, state_col)
    
    treatment_col = df.columns.map(treatment_mapper)
    treatment_list = list(treatment_col.unique())
    treatment_code = map_values_to_value_list(treatment_list, treatment_col)

    midx = pd.MultiIndex(levels=[cell_line_list, state_list, treatment_list], 
                     codes=[cell_line_code, state_code, treatment_code],
                     names=["cell_line", "state", "treatment"])
    return midx


def diffacto_to_midx_df(diffacto, protein_index = True):
    if protein_index == True:
        diffacto_vals = diffacto.iloc[:,diffacto.columns.isin(['Protein','RKO_D_0', 'RKO_D_1',
               'RKO_D_2', 'RKO_D_3', 'RKO_D_4', 'RKO_D_5', 'RKO_D_6', 'RKO_D_7',
               'RKO_D_8', 'RKO_D_9', 'RKO_S_0', 'RKO_S_1', 'RKO_S_2', 'RKO_S_3',
               'RKO_S_4', 'RKO_S_5', 'RKO_S_6', 'RKO_S_7', 'RKO_S_8', 'RKO_S_9',
               'A549_D_0', 'A549_D_1', 'A549_D_2', 'A549_D_3', 'A549_D_4', 'A549_D_5',
               'A549_D_6', 'A549_D_7', 'A549_D_8', 'A549_D_9', 'A549_S_0', 'A549_S_1',
               'A549_S_2', 'A549_S_3', 'A549_S_4', 'A549_S_5', 'A549_S_6', 'A549_S_7',
               'A549_S_8', 'A549_S_9', 'MCF7_D_0', 'MCF7_D_1', 'MCF7_D_2', 'MCF7_D_3',
               'MCF7_D_4', 'MCF7_D_5', 'MCF7_D_6', 'MCF7_D_7', 'MCF7_D_8', 'MCF7_D_9',
               'MCF7_S_0', 'MCF7_S_1', 'MCF7_S_2', 'MCF7_S_3', 'MCF7_S_4', 'MCF7_S_5',
               'MCF7_S_6', 'MCF7_S_7', 'MCF7_S_8', 'MCF7_S_9'])].set_index("Protein")
    else:
        diffacto_vals = diffacto.iloc[:,diffacto.columns.isin(['RKO_D_0', 'RKO_D_1',
       'RKO_D_2', 'RKO_D_3', 'RKO_D_4', 'RKO_D_5', 'RKO_D_6', 'RKO_D_7',
       'RKO_D_8', 'RKO_D_9', 'RKO_S_0', 'RKO_S_1', 'RKO_S_2', 'RKO_S_3',
       'RKO_S_4', 'RKO_S_5', 'RKO_S_6', 'RKO_S_7', 'RKO_S_8', 'RKO_S_9',
       'A549_D_0', 'A549_D_1', 'A549_D_2', 'A549_D_3', 'A549_D_4', 'A549_D_5',
       'A549_D_6', 'A549_D_7', 'A549_D_8', 'A549_D_9', 'A549_S_0', 'A549_S_1',
       'A549_S_2', 'A549_S_3', 'A549_S_4', 'A549_S_5', 'A549_S_6', 'A549_S_7',
       'A549_S_8', 'A549_S_9', 'MCF7_D_0', 'MCF7_D_1', 'MCF7_D_2', 'MCF7_D_3',
       'MCF7_D_4', 'MCF7_D_5', 'MCF7_D_6', 'MCF7_D_7', 'MCF7_D_8', 'MCF7_D_9',
       'MCF7_S_0', 'MCF7_S_1', 'MCF7_S_2', 'MCF7_S_3', 'MCF7_S_4', 'MCF7_S_5',
       'MCF7_S_6', 'MCF7_S_7', 'MCF7_S_8', 'MCF7_S_9'])]


    midx = diffacto_col_to_mIdx(diffacto_vals)
    diffacto_vals.columns = midx
    return diffacto_vals
