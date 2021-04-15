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

import scipy
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


def drop_zero_row(df):
    return df[(df.T != 0).any()]
    
# filter df
df_base.columns
df = df_raw[df_raw.PEP < 0.05] # 5% PEP removed 9578 peptides
df = df[df["Missed cleavages"] < 3] # 0 removed
df = df.set_index("Leading razor protein")

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

#kde_matrix_plot_batch(df_irs_tmm, log = False, suptitle = "Protein(Raw-Aitchison-Top3-SL-IRS-TMM) (Channel intensities together)", legend_size = 6)
kde_matrix_plot_batch(df_protQuant, log = False, suptitle = "Protein(Raw-SL-IRS-TMM-Top3) (Channel intensities together)", legend_size = 6)


kde_matrix_plot_batch(df_irs_tmm*df_signal, log = False, suptitle = "raw -> aitchison -> SL -> IRS -> TMM (Channel intensities together)", legend_size = 6)
kde_matrix_plot_batch(data_aitchison, log = False, suptitle = "raw -> SL -> IRS -> TMM -> Aitchison (Channel intensities together)", legend_size = 6)


##############
# AITCHISOLN #
##############

data_aitchison = aitchison_transform(df_int.replace(0,np.nan)) #this is better normalized... plot the kde for this in report

data_aitchison = aitchison_transform(df_irs_tmm)

kde_matrix_plot_all_channels(data_aitchison, log = False, suptitle ="Aitchison")
kde_matrix_plot_batch(data_aitchison, log = False, suptitle = "Aitchison (Channel intensities together)", legend_size = 6)
      


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
    df_pca = df_int.iloc[:,df_int.columns.get_level_values(seperator) == seperate_by]

    df_pca = df_pca.fillna(0).T # transpose because other did not work
    
    features = df_int.index.values
    
    # Separating out the features
    x = df_pca.loc[:, features].values
    
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
    plt.setp(axes[row,col].get_legend().get_texts(), fontsize='6')


fig, axes = plt.subplots(3, 3, figsize = (30,24))
cell_lines = ["A549", "RKO", "MCF7"]
color_markers = [["state", "replicate"], ["treatment", "state"], ["treatment", "replicate"]]
for row in range(3):
    color_marker = color_markers[row]
    for col in range(3):
        cell_line = cell_lines[col]
        pca_plot_ax(axes, row, col, data_aitchison, classification = color_marker[0], marker = color_marker[1], seperator = "cell_line",
                    seperate_by = cell_line, title = None)
plt.suptitle("Raw->SL->IRS->TMM->Aitchison")

fig, axes = plt.subplots(3, 3, figsize = (30,24))
cell_lines = ["A549", "RKO", "MCF7"]
color_markers = [["state", "replicate"], ["treatment", "state"], ["treatment", "replicate"]]
for row in range(3):
    color_marker = color_markers[row]
    for col in range(3):
        cell_line = cell_lines[col]
        pca_plot_ax(axes, row, col, df_irs_tmm, classification = color_marker[0], marker = color_marker[1], seperator = "cell_line",
                    seperate_by = cell_line, title = None)
plt.suptitle("Raw->SL-IRS->TMM")



#########################
# Protein summarization #
#########################
        
df = data_aitchison

df = df_int
df = pd.DataFrame(df.values, index = df.index, columns = data_aitchison.columns.get_level_values("experiment"))

protein = "P55011"

df[df.index == protein]

import time
start = time.time()
#protein_array = []- same a df.index.unique()
protein_quants = []
experiment_array = []
for col_i in range(len(df.columns)):
    protein_quant_array = []
    for protein in df.index.unique():
        top3 = df[df.index == protein].iloc[:,col_i].nlargest(3)
        if len(top3)>2:
            proteinQuant = top3.mean()
        else:
            proteinQuant = np.nan
        protein_quant_array.append(proteinQuant)
    print(time.time()-start)
    protein_quants.append(protein_quant_array)
    experiment_array.append(top3.name)

df_protQuant = pd.DataFrame(protein_quants, index = experiment_array, columns = df.index.unique()).T
df_protQuant.to_csv("PROT_QUANT_DF_aitchison_unnormalized.csv", sep = "\t")
end = time.time()
    
print(end-start)
df_protQuant = pd.read_csv("PROT_QUANT_DF.csv", sep = "\t", index_col = 0)
df_protQuant = pd.read_csv("PROT_QUANT_DF.csv", sep="\t", index_col=0, header=[0,1,2,3,4,5,6])

nlargest(3, columns = "0_A549_D_Rep1")
df[df.index == protein].nlargest(3, columns = "1_A549_D_Rep1" )
df[df.index == protein].nlargest(3, columns = "9_RKO_S_Rep3")
df[df.index == protein].nlargest(3, columns = df.columns)


###########################
# Protein sum col to mIdx #
###########################

def protSum_col_to_mIdx(df_int):
    treatment_mapper = lambda x : x.split("_")[0]
    specie_mapper = lambda x : x.split("_")[1]
    state_mapper = lambda x : x.split("_")[2]
    replicate_mapper = lambda x : x.split("_")[3]
    batch_mapper = lambda x : x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[3]
    experiment_mapper = lambda x : x.split("_")[0] + "_" + x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[3]
    sample_mapper = lambda x : x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[0]

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



def protSum_intensities_to_midx_df(df_int):
    midx = protSum_col_to_mIdx(df_int)
    df_int.columns = midx
    return df_int

df_protQuant = protSum_intensities_to_midx_df(df_protQuant)

##################
# Plot protQuant #
##################

# Do I need to normalize again?
# Do I need to 


kde_matrix_plot_all_channels(df_protQuant, log=False, suptitle = "Protein quantities")
kde_matrix_plot_batch(df_protQuant, log=False, suptitle = "Protein quantities") #looks ok?


######################
# Plot PCA protQuant #
######################

fig, axes = plt.subplots(3, 3, figsize = (30,24))
cell_lines = ["A549", "RKO", "MCF7"]
color_markers = [["state", "replicate"], ["treatment", "state"], ["treatment", "replicate"]]
for row in range(3):
    color_marker = color_markers[row]
    for col in range(3):
        cell_line = cell_lines[col]
        pca_plot_ax(axes, row, col, df_protQuant, classification = color_marker[0], marker = color_marker[1], seperator = "cell_line",
                    seperate_by = cell_line, title = None)
plt.suptitle("Raw->SL-IRS->TMM->Aitchison->Top3")

# This looks, very good?

###################################
# Protein average between samples #
###################################

df_protQuant = pd.read_csv("PROT_QUANT_DF.csv", sep = "\t", index_col = 0)
df_protQuant = protSum_intensities_to_midx_df(df_protQuant)

def aggregate_protein_quantity(df_protQuant):
    """
    Get the mean aggregate of protein quantites.
    """
    df_quant = pd.DataFrame()
    for sample in df_protQuant.columns.get_level_values("sample").unique():
        df_sample = df_protQuant.iloc[:,df_protQuant.columns.get_level_values("sample") == sample]
        #count = df_sample.count(axis=1) > 1 #sample treshold
        #df_sample = (df_sample.T*count.values).T
        df_sample_quant = pd.DataFrame(df_sample.mean(axis=1), columns = [sample])
        df_quant = pd.concat([df_quant, df_sample_quant], axis = 1)
    return df_quant.replace(0, np.nan)

def get_p_matrix(df_protQuant):
    df_pVals = pd.DataFrame()
    df_tStats = pd.DataFrame()
    for treatment in df_protQuant.columns.get_level_values("treatment").unique()[1:]:
        df_control = df_protQuant.iloc[:,df_protQuant.columns.get_level_values("treatment") == "0"]    
        for cell_line in df_protQuant.columns.get_level_values("cell_line").unique():
            for state in df_protQuant.columns.get_level_values("state").unique():
                sample = cell_line + "_" + state + "_" + treatment
                df_sample = df_protQuant.iloc[:, df_protQuant.columns.get_level_values("sample") == sample]
                df_sample_control = df_control[cell_line][state]
                tStats = pd.DataFrame(scipy.stats.ttest_ind(df_sample,df_sample_control, axis = 1)[0], columns = [sample], index = df_sample.index)
                pVals = pd.DataFrame(scipy.stats.ttest_ind(df_sample,df_sample_control, axis = 1)[1], columns = [sample], index = df_sample.index)
                df_tStats = pd.concat([df_tStats, tStats], axis = 1)
                df_pVals = pd.concat([df_pVals, pVals], axis = 1)
    return df_pVals, df_tStats
    
    
df_pVals, df_tStats = get_p_matrix(df_protQuant)

df_quant = aggregate_protein_quantity(df_protQuant)

df_quant = diffacto_to_midx_df(df_quant, protein_index=False) # Function is below at Diffacto section



plot_diffacto_pca(df_quant) # JUMP DOWN TO DIFFACTO SECTION TO GET PLOTTING FUNCTION....
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (all)")
plot_diffacto_pca_cell_line(df_quant.A549) # Good !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (A549)")
plot_diffacto_pca_cell_line(df_quant.RKO) # Good  !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (RKO)")
plot_diffacto_pca_cell_line(df_quant.MCF7) # Good !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (MCF7)")


#######################################################
# DIFFERENCE BETWEEN RESPONCS OF ALIVE AND DEAD CELLS #
#######################################################


def get_log2FC_regulation(df_quant):
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

df_log2FC, df_pVals, df_tStats = get_log2FC(df_protQuant)


################
# VOLCANO PLOT #
################



def get_significant_proteins(df_log2FC, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05):
    df_volcano = pd.DataFrame([df_log2FC[cell_line][state][treatment], df_pVals[cell_line][state][treatment]], index = ["log2FC", "p-value"]).T
    df_volcano = df_volcano.dropna()
    
    def volcano_hue(log2fc, pVal, fc_treshold = 1.0, pVal_treshold = 0.05):
        if log2fc > fc_treshold:
            if -np.log10(pVal) > -np.log10(pVal_treshold):
                return "up-regulated"
            else:
                return "non-significant"
        if log2fc < -fc_treshold:
            if -np.log10(pVal) > -np.log10(pVal_treshold):
                return "down-regulated"
            else:
                return "non-significant"
        else:
            return "non-significant"    
    
    #fc_treshold = 1.0
    #pVal_trehsold = 0.05
    
    df_volcano["-log10(p-value)"] = -np.log10(df_volcano["p-value"])
    df_volcano["sign"] = df_volcano.apply(lambda x: volcano_hue(x["log2FC"], x["p-value"], fc_treshold = fc_treshold, pVal_treshold = pVal_treshold), axis = 1)
    return df_volcano

def volcano_plot(df_log2FC, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05):
    df_volcano = get_significant_proteins(df_log2FC, cell_line, state, treatment, fc_treshold = fc_treshold, pVal_treshold = pVal_treshold)
    sns.scatterplot(x="log2FC", y="-log10(p-value)", hue = "sign", data = df_volcano)
    for i in range(np.shape(df_volcano[df_volcano.sign != "non-significant"])[0]):
        vals = df_volcano[df_volcano.sign != "non-significant"].iloc[i,:]
        plt.text(vals["log2FC"], -np.log10(vals["p-value"]), vals.name , 
                 horizontalalignment='left', size='small', color='black', weight='light')
    plt.show()


cell_line = "A549"
state = "D"
treatment = "1"
df_significant = get_significant_proteins(df_log2FC, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05)
volcano_plot(df_log2FC, cell_line, state, treatment, fc_treshold = 0.5, pVal_treshold = 0.05)



##########################
# UP AND DOWN REGULATION #
##########################


df_regulation = get_log2FC_regulation(df_quant)


sign_count_treshold = 14 # we want more than half to be regualted in the same way.


upRegulation=True

def get_mean_regulation_per_state(df_quant, state, up_regulation, fc_treshold, sign_count_treshold = 14):
    df_regulation = get_log2FC_regulation(df_quant)
    df_res = df_regulation.iloc[:,df_regulation.columns.get_level_values("state") == state]
    
    if up_regulation == True:
        # up-regulation
        treshold_proteins = ((df_res > 0).sum(axis = 1) >sign_count_treshold)
        df_res = df_res[treshold_proteins].mean(axis = 1)
        df_res = df_res.sort_values(ascending=False)
        df_res = df_res[df_res > fc_treshold] #more than fc_treshold on up-regulation
    else:
        # down-regulation
        treshold_proteins = ((df_res < 0).sum(axis = 1) >sign_count_treshold)
        df_res = df_res[treshold_proteins].mean(axis = 1)
        df_res = df_res.sort_values(ascending=True)
        df_res = df_res[df_res < fc_treshold] #less than fc_treshold on down-regulation
    return df_res    

treshold = 0.2
df_up_S = get_mean_regulation_per_state(df_quant, state = "S", up_regulation = True, fc_treshold = treshold, sign_count_treshold = 14)
df_down_S = get_mean_regulation_per_state(df_quant, state = "S", up_regulation = False, fc_treshold = -treshold, sign_count_treshold = 14) 
df_up_D = get_mean_regulation_per_state(df_quant, state = "D", up_regulation = True, fc_treshold = treshold,sign_count_treshold = 14)
df_down_D = get_mean_regulation_per_state(df_quant, state = "D", up_regulation = False, fc_treshold = -treshold, sign_count_treshold = 14)             

# Make venn diagram of these...
from matplotlib_venn import venn2

venn2([set(df_up_S.index), set(df_up_D.index)], ("Surviving cells", "Dead cells"))
plt.title("Up-regulated, FC = " + str(treshold))
df_up_S.head(20).plot.bar(title = "up-regulated surviving, FC = 0.2, top 20")
df_up_D.head(20).plot.bar(title = "up-regulated dead, FC = 0.2, top 20")

venn2([set(df_down_S.index), set(df_down_D.index)], ("Surviving cells", "Dead cells"))
plt.title("Down-regulated, FC = " + str(treshold))
df_down_S.head(20).plot.bar(title = "down-regulated surviving, FC = -0.2, top 20")
df_down_D.head(20).plot.bar(title = "down-regulated dead, FC = -0.2, top 20")



##################################
# VENN with significant proteins #
##################################
df_log2FC, df_pVals, df_tStats = get_log2FC(df_protQuant)
cell_line = "RKO"
treatment = "2"
fc_tresh = 1.0
pVal_tresh = 0.05


for treatment in range(1,10):
    i = treatment
    treatment = str(treatment)
    plt.subplot(3,3,i)
    df_significant_S = get_significant_proteins(df_log2FC, df_pVals, cell_line, "S", treatment, fc_treshold = fc_tresh, pVal_treshold = pVal_tresh)
    df_significant_D = get_significant_proteins(df_log2FC, df_pVals, cell_line, "D", treatment, fc_treshold = fc_tresh, pVal_treshold = pVal_tresh)
    
    df_S = df_significant_S[df_significant_S.sign != "non-significant"]
    df_D = df_significant_D[df_significant_D.sign != "non-significant"]
    
    venn2([set(df_S.index), set(df_D.index)], ("Surviving cells", "Dead cells"))
    plt.title("treatment " + treatment)
plt.suptitle("Significant proteins, " + cell_line + ", FC = " + str(fc_tresh) + ", p-value = " + str(pVal_tresh))
############
#
#############

 
 
######################
######################
# DIFFACTO ###########
######################
######################


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



###########
# PCA #####
###########

# For a more detailed PCA check ws_202104066_withPCA.py

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import seaborn as sns


def plot_diffacto_pca(df):
    marker = "state"
    colour = "cell_line"
    
    df_pca = df.T
    df_pca = df_pca.fillna(0) # Check if this is needed.
    features = df.index.values
    x = df_pca.loc[:, features].values
    x = StandardScaler().fit_transform(x)
    
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ["principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                              'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1])])
    finalDf = principalDf
    finalDf[colour] = df_pca.index.get_level_values(colour)
    finalDf[marker] = df_pca.index.get_level_values(marker)
    
    sns.scatterplot(data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                    y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]), 
                    hue = colour, style = marker, s = 100)

def plot_diffacto_pca_cell_line(df):
    marker = "state"
    #colour = "cell_line"
    
    df_pca = df.T
    df_pca = df_pca.fillna(0) # Check if this is needed.
    features = df.index.values
    x = df_pca.loc[:, features].values
    x = StandardScaler().fit_transform(x)
    
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ["principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                              'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1])])
    finalDf = principalDf
    #finalDf[colour] = df_pca.index.get_level_values(colour)
    finalDf[marker] = df_pca.index.get_level_values(marker)
    
    sns.scatterplot(data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                    y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]),
                    style = marker, s = 100)
df = diffacto_to_midx_df(diffacto)

plot_diffacto_pca(df)
plt.title("PCA (Aitchison -> Diffacto)")


