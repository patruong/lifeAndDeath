#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 00:30:17 2021

@author: ptruong
"""

import os 
import pandas as pd
import numpy as np

from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt


from get_columns import get_cell_line_state_replicate, get_base_cols_peptide, get_all_reporter_intensity_correct, get_reporter_intensity_without_control
from column_mapper import col_to_treatment_mapper, treatment_nomenclature_map_dict, col_to_cell_line_mapper, col_to_state_mapper, col_to_rep_mapper
pd.options.mode.chained_assignment = None
os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")
df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")


cell_lines, states, replicates = get_cell_line_state_replicate()
base_cols = get_base_cols_peptide()
reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()
df_base = df[get_base_cols_peptide()]



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
    treatment_mapper = lambda x : x.split(" ")[3]
    specie_mapper = lambda x : x.split(" ")[4].split("_")[0]
    state_mapper = lambda x : x.split(" ")[4].split("_")[1]
    replicate_mapper = lambda x : x.split(" ")[4].split("_")[2]
    
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
    
    midx = pd.MultiIndex(levels=[specie_list, state_list, treatment_list, replicate_list], 
                         codes=[specie_code, state_code, treatment_code, replicate_code])
    df_int.columns = midx
    return df_int

# filter df
df_base.columns
df = df_raw[df_raw.PEP < 0.05] # 5% PEP removed 9578 peptides
df = df[df["Missed cleavages"] < 3] # 0 removed

df_int = df[reporter_intensity_corrected_cols]
df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
df_int = intensities_to_midx_df(df_int)

df_int.T





