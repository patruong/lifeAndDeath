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
###########################################
# RAW -- >  NORM --> COMBAT --> AITCHISON #
###########################################



########
# NORM #
########
    
def norm_SL(df_int):
    target = df_int.sum().mean()
    norm_fac = target/df_int.sum()
    df_norm = df_int*norm_fac
    return df_norm

df_norm = norm_SL(df_int)
kde_matrix_plot_all_channels(df_norm, log = False, suptitle ="log2(norm)")
kde_matrix_plot_batch(df_norm, log = False, suptitle = "log2(Norm intensity) (Channel intensities together)", legend_size = 6)


#######
# TMM #
#######

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

sl_tmm = pd.Series(sl_tmm[0], df_norm.columns)

df_tmm = df_norm/(sl_tmm)

kde_matrix_plot_all_channels(df_tmm, log = True, suptitle ="log2(TMM)")
kde_matrix_plot_batch(df_tmm, log = True, suptitle = "log2(TMM) (Channel intensities together)", legend_size = 6)

#######
# IRS #
#######

irs = pd.DataFrame()
for batch in df_norm.columns.get_level_values("batch").unique():
    irs_sum = df_norm.iloc[:,df_norm.columns.get_level_values("batch") == batch].sum(axis = 1)
    irs_sum = pd.DataFrame(irs_sum.values, columns = ["sum_" + batch], index = df_norm.index)
    irs = pd.concat([irs, irs_sum], axis = 1)

irs_geoAvg = np.exp(np.log(df_norm.replace(0, np.nan)).mean(axis=1))
irs_geoAvg = pd.DataFrame(irs_geoAvg, columns = ["geoAvg"])
irs = pd.concat([irs, irs_geoAvg], axis = 1)
irs = irs.replace(0,np.nan)

i = 0
for batch in df_norm.columns.get_level_values("batch").unique():
    irs_fac = irs.geoAvg/irs.iloc[:,i]
    irs_fac = pd.DataFrame(irs_fac.values, columns = ["fac_"+batch], index = df_norm.index)
    irs = pd.concat([irs, irs_fac], axis = 1)    
    i+=1

df_irs = pd.DataFrame()
for batch in df_norm.columns.get_level_values("batch").unique():
    exp_sl = df_norm.iloc[:,df_norm.columns.get_level_values("batch") == batch]
    irs_fac = pd.DataFrame(irs["fac_" + batch])
    
    data_irs = irs_fac.values*exp_sl
    df_irs = pd.concat([df_irs, data_irs], axis = 1)


kde_matrix_plot_all_channels(df_irs, log = True, suptitle ="log2(IRS)")
kde_matrix_plot_batch(df_irs, log = True, suptitle = "log2(IRS) (Channel intensities together)", legend_size = 6)

#############
# IRS - TMM #
#############


irs_tmm = calcNormFactors(df_irs.fillna(0))

irs_tmm = pd.Series(irs_tmm[0], df_irs.columns)

df_irs_tmm = df_irs/(irs_tmm)


kde_matrix_plot_all_channels(df_irs_tmm, log = True, suptitle ="log2(IRS-TMM)")
kde_matrix_plot_batch(df_irs_tmm, log = True, suptitle = "log2(IRS-TMM) (Channel intensities together)", legend_size = 6)



##############
# AITCHISOLN #
##############

data_aitchison = aitchison_transform(df_irs_tmm)

kde_matrix_plot_all_channels(data_aitchison, log = False, suptitle ="Aitchison")
kde_matrix_plot_batch(data_aitchison, log = False, suptitle = "Aitchison (Channel intensities together)", legend_size = 6)




#################################
# Generate diffacto input files #
#################################


data_aitchison["sequence"] = df_raw.Sequence
data_aitchison = data_aitchison.set_index("sequence")

df_signal["sequence"] = df_raw.Sequence
df_signal = df_signal.set_index("sequence")

data_aitchison_true_signal = data_aitchison*df_signal

diffacto_input = pd.DataFrame(data_aitchison_true_signal.values, index = data_aitchison_true_signal.index.values, columns = data_aitchison_true_signal.columns.get_level_values("experiment").values)


np.log2(diffacto_input).to_csv("peptide_tryptic_aitchison_diffacto_input.csv", sep = ",")

sample_rep_name = data_aitchison_true_signal.columns.get_level_values("experiment").values
sample_name = data_aitchison_true_signal.columns.get_level_values("sample").values

pd.DataFrame([sample_rep_name, sample_name]).T.to_csv("peptide_tryptic_aitchison_diffacto_input_sample.lst", index = False, header = False, sep = "\t")
diffacto_input = pd.read_csv("peptide_tryptic_aitchison_diffacto_input.csv", sep = ",", index_col = 0)

def check_if_row_elements_equal(df):
    return df.eq(df.iloc[:, 0], axis=0)


def to_diffacto_input(data_aitchison, df_signal):

    data_aitchison["sequence"] = df_raw.Sequence
    
    data_aitchison = data_aitchison.set_index("sequence")

    df_signal["sequence"] = df_raw.Sequence

    df_signal = df_signal.set_index("sequence")
    #print(df_signal)
    data_aitchison_true_signal = data_aitchison*df_signal
    
    diffacto_input = pd.DataFrame(data_aitchison_true_signal.values, index = data_aitchison_true_signal.index.values, columns = data_aitchison_true_signal.columns.get_level_values("experiment").values)
    
    sample_rep_name = data_aitchison_true_signal.columns.get_level_values("experiment").values
    sample_name = data_aitchison_true_signal.columns.get_level_values("sample").values
    diffacto_sample_list = pd.DataFrame([sample_rep_name, sample_name]).T
    return diffacto_input, diffacto_sample_list
    
diffacto_input, diffacto_sample_list = to_diffacto_input(data_aitchison, df_signal)
diffacto_input.to_csv("peptide_tryptic_aitchison_diffacto_input.csv", sep = ",")
diffacto_sample_list.to_csv("peptide_tryptic_aitchison_diffacto_input_sample.lst", index = False, header = False, sep = "\t")                        

pd.read_csv("peptide_tryptic_aitchison_diffacto_input.csv")
#####################
# read in diffacto ##
#####################
    
diffacto = pd.read_csv("diffacto_output.protein.txt", sep = "\t")
diffacto_FDR = pd.read_csv("diffacto_output.protein.FDR", sep = "\t")
diffacto["MCFDR"] = diffacto_FDR["MCFDR"]
diffacto = diffacto[diffacto["MCFDR"] < 0.05] #treshold on monte-carlo FDR
diffacto = diffacto[diffacto["S/N"] > 0] # treshold on signal-to-noise ratio


############
# DIFFACTO #
############

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


def diffacto_to_midx_df(diffacto):
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

    midx = diffacto_col_to_mIdx(diffacto_vals)
    diffacto_vals.columns = midx
    return diffacto_vals


df = diffacto_to_midx_df(diffacto)




###########
# PCA #####
###########

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import seaborn as sns



def pca_plot_ax(axes, row, col, df_int, classification = "state", marker = "replicate", seperator = "cell_line", seperate_by = "A549",
            title = None):
    # splice df_int for pca now
    #df_pca = df_int[df_int[seperator] == seperate_by]
    df_pca = df_int.iloc[:,df_int.columns.get_level_values(seperator) == seperate_by]
    #df_pca = df_pca.fillna(0)
    df_pca = df_pca.fillna(0).T # transpose because other did not work
    #weights = df_pca.notna()*1
    
    features = df_int.index.values
    
    # Separating out the features
    x = df_pca.loc[:, features].values
    
    # Separating out the target
    # classification = "state # Choose target - treatment, cell_line, state here
    # marker = "rep" #marker
    # y = df_pca.loc[:,[classification]].values
    
    # Standardizing the features
    x = StandardScaler().fit_transform(x)
    
    # Missing value impuration
    imputer = SimpleImputer(missing_values=np.nan, strategy="constant", fill_value = 0)
    x = imputer.fit_transform(x)
    
    # pca
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ["principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                              'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1])])
    #finalDf = pd.concat([principalDf, df_pca[[classification]].reset_index()[classification]], axis = 1)
    finalDf = principalDf
    finalDf[classification] = df_pca.index.get_level_values(classification)
    if marker != False:
        #finalDf = pd.concat([finalDf, df_pca[[marker]].reset_index()[marker]], axis = 1)
        finalDf[marker] = df_pca.index.get_level_values(marker)
    #Visualization
    #fig, ax = plt.subplots(figsize=(16,10))
    if marker != False:
        sns.scatterplot(ax=axes[row,col] , data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                        y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]), 
                        hue = classification, style = marker, s = 100)
    else:
        sns.scatterplot(ax=axes[row,col], data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]), 
                        y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]), 
                hue = classification, s = 100)
    axes[row, col].set_title(seperate_by)


fig, axes = plt.subplots(3, 3, figsize = (30,24))

cell_lines = ["A549", "RKO", "MCF7"]
color_markers = [["state", "replicate"], ["treatment", "state"], ["treatment", "replicate"]]
for row in range(3):
    color_marker = color_markers[row]
    for col in range(3):
        cell_line = cell_lines[col]
        pca_plot_ax(axes, row, col, df_norm, classification = color_marker[0], marker = color_marker[1], seperator = "cell_line",
                    seperate_by = cell_line, title = None)





