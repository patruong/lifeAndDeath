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
df_base = df_raw[get_base_cols_peptide()]


def col_to_mIdx(df_int):
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

def drop_zero_row(df):
    return df[(df.T != 0).any()]
    
# filter df
df_base.columns
df = df_raw[df_raw.PEP < 0.05] # 5% PEP removed 9578 peptides
df = df[df["Missed cleavages"] < 3] # 0 removed

df_int = df[reporter_intensity_corrected_cols]
df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
df_int = drop_zero_row(df_int) # dropped last zero-row
midx = col_to_mIdx(df_int)

# Raw
df_int = intensities_to_midx_df(df_int)

# log2
df_int = df_int.replace(0, np.nan)
df_int = np.log2(df_int)


# Aitchison
#df_int = closure(df_int)
#df_int = pd.DataFrame(df_int, columns = midx)

# Aitchison multiplicative_replacement
#df_int = multiplicative_replacement(df_int)
#df_int = pd.DataFrame(df_int, columns = midx)

#######################
# BATCH NORMALIZATION #
# Multiplication is not the same in log-space so we need to normalize before any transformation....
#######################

df_int.sum() #check intensity sum of each channel

target = df_int.sum().mean()
norm_fac = target/df_int.sum()
df_norm = df_int*norm_fac

df_norm.sum() #check again intensity sum of each channel.

df_norm_log2 = df_norm.replace(0, np.nan)
df_norm_log2 = np.log2(df_norm_log2)


#import statsmodels.robust.norms.TrimmedMean
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def calcNormFactors(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
      r_from_pd_df = ro.conversion.py2rpy(df)
    
    edgeR = importr("edgeR")
    raw_tmm = edgeR.calcNormFactors(r_from_pd_df)
    
    with localconverter(ro.default_converter + pandas2ri.converter):
      pd_from_r_df = ro.conversion.rpy2py(raw_tmm)
      
    raw_tmm_pd = pd.DataFrame(pd_from_r_df, index =midx).T.values
    return raw_tmm_pd

sl_tmm = calcNormFactors(df_norm)

df_tmm = df_norm/sl_tmm
df_tmm = df_tmm.replace(0, np.nan)
df_tmm = np.log2(df_tmm)

# IRS

df_geoMean  = np.exp(np.log(df_norm.replace(0,np.nan)).mean(axis=1))
irs_fac = df_geoMean / df_norm.sum(axis=1) 
df_irs =(df_norm.T*irs_fac).T

irs_tmm = calcNormFactors(df_irs)
df_irs_tmm = df_irs/irs_tmm

a = pd.DataFrame([[2,3,4]]).T
b = pd.DataFrame([[1,2,3]]).T

a = pd.DataFrame([[1,2,3],[2,3,4],[5,6,7]])
b = pd.DataFrame([[1,2,3]]).values
########
# PLOT #
########

def plot_intensity_histogram(df_int, min_x = 0, max_x = 30, step = 0.1, title = "title"):
    """
    df_int with midx
    
    plots the step histogram for {cell_line}_{state}
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    col_i = 0
    for cell_line in cell_lines:
        for state in states:
            df_state = df_int.iloc[:, df_int.columns.get_level_values(level=1) == state]
            df_cell_state = df_state.iloc[:, df_state.columns.get_level_values(level=0) == cell_line]
            plot_label = True
            for i in range(np.shape(df_cell_state)[1]):
                if plot_label == True:
                    plt.hist(df_cell_state.iloc[:,i] , bins=np.arange(min_x,max_x, step), 
                             histtype="step", color = colors[col_i], alpha = 0.4,
                             label = cell_line + "_" + state)
                    plot_label = False
                else:
                    plt.hist(df_cell_state.iloc[:,i] , bins=np.arange(min_x,max_x, step), 
                             histtype="step", color = colors[col_i], alpha = 0.4)
            col_i+=1
    plt.title(title)
    plt.legend()

plot_intensity_histogram(np.log2(df_int.replace(0,np.nan)), min_x = 0, max_x = 25, step = 0.1, title = "raw")
plot_intensity_histogram(np.log2(df_norm.replace(0,np.nan)), min_x = 0, max_x = 25, step = 0.1, title = "Batch normalization - equal signal per channel")
plot_intensity_histogram(df_tmm, min_x = 0, max_x = 25, step = 0.1, title = "TMM")
plot_intensity_histogram(np.log2(df_irs.replace(0,np.nan)), min_x = 0, max_x = 25, step = 0.1, title = "IRS")
plot_intensity_histogram(np.log2(df_irs_tmm.replace(0,np.nan)), min_x = 0, max_x = 25, step = 0.1, title = "IRS-TMM")

def plot_intensity_boxplot(df_int, title = "title"):
    """
    Boxplot details.
    # https://stackoverflow.com/questions/41997493/python-matplotlib-boxplot-color

    df_int with midx
    
    plots the step histogram for {cell_line}_{state}
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    col_i = 0
    pos_i = 0
    for cell_line in cell_lines:
        for state in states:
            df_state = df_int.iloc[:, df_int.columns.get_level_values(level=1) == state]
            df_cell_state = df_state.iloc[:, df_state.columns.get_level_values(level=0) == cell_line]
            plot_label = True
            for i in range(np.shape(df_cell_state)[1]):
                if plot_label == True:
                    plt.boxplot(df_cell_state.iloc[:,i].dropna(), positions = [pos_i], patch_artist=True, 
                                notch=True, medianprops=dict(color="grey"),
                                boxprops=dict(facecolor = colors[col_i]),)
                                # No label yet, so if-else does not matter.......
                    plot_label = False
                else:
                    plt.boxplot(df_cell_state.iloc[:,i].dropna(), positions = [pos_i], patch_artist=True, 
                                notch=True, medianprops=dict(color="grey"),
                                boxprops=dict(facecolor = colors[col_i]))
                pos_i+=1
            col_i+=1
    plt.title(title)
    plt.legend()

plot_intensity_boxplot(np.log2(df_int.replace(0,np.nan)), title = "raw")
plot_intensity_boxplot(np.log2(df_norm.replace(0,np.nan)), title = "Batch normalization - equal signal per channel")
plot_intensity_boxplot(df_tmm, title = "TMM")
plot_intensity_boxplot(np.log2(df_irs.replace(0,np.nan)), title = "IRS")
plot_intensity_boxplot(np.log2(df_irs_tmm.replace(0,np.nan)), title = "IRS-TMM")

