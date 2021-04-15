#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 01:58:31 2021

@author: ptruong
"""


import os 
import time

import pandas as pd
import numpy as np

from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

from combat.pycombat import pycombat

import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from get_columns import get_cell_line_state_replicate, get_base_cols_peptide, get_all_reporter_intensity_correct, get_reporter_intensity_without_control
from column_mapper import col_to_treatment_mapper, treatment_nomenclature_map_dict, col_to_cell_line_mapper, col_to_state_mapper, col_to_rep_mapper
pd.options.mode.chained_assignment = None
os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")
df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")


#cell_lines, states, replicates = get_cell_line_state_replicate()
base_cols = get_base_cols_peptide()
reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()

cell_lines, states, replicates = [], [], []
for x in reporter_intensity_corrected_cols:
    cell_line = x.split(" ")[-1].split("_")[0]
    state = x.split(" ")[-1].split("_")[1]
    replicate = x.split(" ")[-1].split("_")[2]
    cell_lines.append(cell_line), states.append(state), replicates.append(replicate)
cell_lines = list(dict.fromkeys(cell_lines).keys())
states = list(dict.fromkeys(states).keys())
replicates = list(dict.fromkeys(replicates).keys())

df_base = df_raw[get_base_cols_peptide()]


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
df_signal = (df_int != 0)


def aitchison_transform_part(df):
    """
    Aitchison tranformation on df with all columns belonging to same batch.
    
    df should consist of all samples tagged together in one channel (i.e. A549_S_rep1 etc.)
    """
    df_aitchison = multiplicative_replacement(df)
    #df_aitchison = closure(df)
    df_idx = df.index
    df_col = df.columns
    df_aitchison = pd.DataFrame(df_aitchison, index = df_idx, columns = df_col)
    return df_aitchison

def aitchison_transform(df_int):
    
    df_aitchison = pd.DataFrame()
    
    start = time.time()
    for cell_line in cell_lines:
        for state in states:
            for rep in replicates:
                df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == ("_".join([cell_line, state, rep]))]
                df_part = df_part[(df_part.T != 0).any()]
                df_part = aitchison_transform_part(df_part)
                df_aitchison = pd.concat([df_aitchison, df_part], axis = 1 )
                print(time.time()-start)
    end=time.time()
    print(end-start)
    return df_aitchison

########
# PLOT #
########

def plot_kde(df_part, log = True, title = "title", legend_size = 6):
    if log == True:
        for i in range(np.shape(df_part)[1]):
            sns.distplot(np.log2((df_part.iloc[:,i]).replace(0,np.nan)), bins = 1000,
                         label = df_part.iloc[:,i].name[-1], kde = True, hist = False)
    else:
        for i in range(np.shape(df_part)[1]):
            sns.distplot((df_part.iloc[:,i]).replace(0,np.nan), bins = 1000,
                         label = df_part.iloc[:,i].name[-1], kde = True, hist = False)
    plt.title(title)
    plt.legend(prop={'size': legend_size})
    
def kde_matrix_plot_all_channels(df_int, log = True, suptitle ="title"):
    i = 0
    for rep in replicates:
        plt.subplot(6, 3, i+1)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["A549", "D", rep]))]
        plot_kde(df_part, log = log, title = "A549_D_"+rep)
        plt.subplot(6, 3, i+2)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["MCF7", "D", rep]))]
        plot_kde(df_part, log = log,  title = "MCF7_D_"+rep)
        plt.subplot(6, 3, i+3)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["RKO", "D", rep]))]
        plot_kde(df_part, log = log,  title = "RKO_D_"+rep)
        
        plt.subplot(6, 3, i+10)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["A549", "S", rep]))]
        plot_kde(df_part, log = log, title = "A549_S_"+rep)
        plt.subplot(6, 3, i+11)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["MCF7", "S", rep]))]
        plot_kde(df_part, log = log, title = "MCF7_S_"+rep)
        plt.subplot(6, 3, i+12)
        df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join(["RKO", "S", rep]))]
        plot_kde(df_part, log = log, title = "RKO_S_"+rep)
        i += 3
    plt.suptitle(suptitle)


def plot_kde_batch(df_int, log = True, title = "title", legend_size = 6):
    if log == True:
        for cell_line in cell_lines:
            for state in states:
                for replicate in replicates:
                    df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join([cell_line, state, replicate]))]
                    sns.distplot(np.log2(df_part.replace(0,np.nan)),
                                 label = cell_line + "_" + state + "_" + replicate,
                                 hist = False, kde = True)
    else:
        for cell_line in cell_lines:
            for state in states:
                for replicate in replicates:
                    df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                              ("_".join([cell_line, state, replicate]))]
                    sns.distplot(df_part,
                                 label = cell_line + "_" + state + "_" + replicate,
                                 hist = False, kde = True)

    plt.title(title)
    plt.legend(prop={'size': legend_size})
    
def kde_matrix_plot_batch(df_int, log = True, suptitle = "title", legend_size = 6):
    plt.subplot(4, 1, 1)
    plot_kde_batch(df_int["A549"], log = log, title = "A549", legend_size = legend_size)
    plt.subplot(4, 1, 2)
    plot_kde_batch(df_int["MCF7"], log = log, title = "MCF7", legend_size = legend_size)
    plt.subplot(4, 1, 3)
    plot_kde_batch(df_int["RKO"], log = log, title = "RKO", legend_size = legend_size)
    plt.subplot(4, 1, 4)
    plot_kde_batch(df_int, log = log, title = "All", legend_size = legend_size)
    plt.suptitle(suptitle)


###########################################
# RAW -- >  NORM --> COMBAT --> AITCHISON #
###########################################

#######
# RAW #
#######
    
kde_matrix_plot_all_channels(df_int, log = True, suptitle ="log2(Raw)")
kde_matrix_plot_batch(df_int, log = True, suptitle = "log2(Raw intensity) (Channel intensities together)", legend_size = 6)



##########
# COMBAT #
##########
    
batch = pd.DataFrame(df_int.columns.get_level_values(4)).T.values
df_combat = pycombat(df_int.fillna(0),batch[0])  #ComBat

kde_matrix_plot_all_channels(df_combat, log = True, suptitle ="log2(ComBat)")
kde_matrix_plot_batch(df_combat*df_signal, log = True, suptitle = "log2(ComBat intensity) (Channel intensities together)", legend_size = 6)



##############
# AITCHISOLN #
##############

data_aitchison = aitchison_transform(df_combat[df_combat>0])

kde_matrix_plot_all_channels(data_aitchison, log = False, suptitle ="Aitchison")
kde_matrix_plot_batch(data_aitchison, log = False, suptitle = "Aitchison (Channel intensities together)", legend_size = 6)


#######################
# AITCHISON ADJUSTED  #
#######################

data_aitchison_adj = (data_aitchison*df_signal).replace(0,np.nan)
data_aitchison_adj = data_aitchison_adj[data_aitchison_adj > 0] # Aitchison > 0, true signal adjusted
                 
kde_matrix_plot_all_channels(data_aitchison_adj, log = False, suptitle ="Aitchison, adjusted")
kde_matrix_plot_batch(data_aitchison_adj, log = False, suptitle = "Aitchison, adjusted (Channel intensities together)", legend_size = 6)

########
# NORM #
########
    
def norm_SL(df_int):
    target = df_int.sum().mean()
    norm_fac = target/df_int.sum()
    df_norm = df_int*norm_fac
    return df_norm

df_norm = norm_SL(data_aitchison_adj)
kde_matrix_plot_all_channels(df_norm, log = False, suptitle ="norm")
kde_matrix_plot_batch(df_norm, log = False, suptitle = "Norm intensity (Channel intensities together)", legend_size = 6)


                        