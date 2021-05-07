#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 12:35:58 2021

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
import scipy.stats as stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

from combat.pycombat import pycombat

import matplotlib.pyplot as plt
import seaborn as sns
import pylab 

os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from get_columns import get_cell_line_state_replicate, get_base_cols_peptide, get_all_reporter_intensity_correct, get_reporter_intensity_without_control
from column_mapper import col_to_treatment_mapper, treatment_nomenclature_map_dict, col_to_cell_line_mapper, col_to_state_mapper, col_to_rep_mapper
from get_variables import get_cell_line_states_replicates_from_reporter_intensity_cols
from midx import col_to_mIdx, intensities_to_midx_df, diffacto_col_to_mIdx, diffacto_to_midx_df
from pd_functions import drop_zero_row
from transform import aitchison_transform_part, aitchison_transform, norm_SL, calcNormFactors, irs_norm
from plot import plot_kde, kde_matrix_plot_all_channels, plot_kde_batch, kde_matrix_plot_batch, plot_intensity_boxplot, plot_diffacto_pca, plot_diffacto_pca_cell_line, pca_plot_ax, get_significant_proteins, volcano_plot, kde_matrix_all_samples 
from top3 import top3, protSum_col_to_mIdx, protSum_intensities_to_midx_df, aggregate_protein_quantity, get_p_matrix
from de_analysis import get_log2FC_regulation, get_log2FC, compute_log2FC, split_surviving_dead
from q_value import qvalues

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.options.mode.chained_assignment = None
os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")
df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")


base_cols = get_base_cols_peptide()

reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()
cell_lines, states, replicates = get_cell_line_states_replicates_from_reporter_intensity_cols(reporter_intensity_corrected_cols)

df_base = df_raw[get_base_cols_peptide()]


# filter df
df = df_raw[df_raw.PEP < 0.05] # 5% PEP removed 9578 peptides
df = df[df["Missed cleavages"] < 3] # 0 removed
df = df.set_index("Leading razor protein")

df.sort_values(by="Score", ascending=False, inplace=True)
df_int = df[reporter_intensity_corrected_cols + ["Score"]]
df_int.sort_values(by="Score", ascending=False, inplace=True)

protein_array = []
proteinQuant_array = []
for protein, peptides in df_int.groupby(df_int.index):
    proteinQuant = peptides.head(1).drop("Score", axis = 1)
    proteinQuant_array.append(proteinQuant.values[0])
    protein_array.append(protein)

df_int = pd.DataFrame(proteinQuant_array, index = protein_array, columns = df_int.drop("Score",axis=1).columns)

df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
df_int = drop_zero_row(df_int) # dropped last zero-row

df_int = intensities_to_midx_df(df_int)


# IRS-TMM NORM
df_norm = norm_SL(df_int)
df_irs = irs_norm(df_norm)
irs_tmm = calcNormFactors(df_irs.fillna(0))
df_irs_tmm = df_irs/(irs_tmm)

df_log2FC = compute_log2FC(df_int)
df_log2FC_norm = compute_log2FC(df_irs_tmm)

S,D,SD = split_surviving_dead(df_log2FC)
Sn,Dn,SDn = split_surviving_dead(df_log2FC_norm)

S,D,SD = split_surviving_dead(df_log2FC, cell_lines = ["RKO"])
plt.figure(figsize=(14,6))
stats.probplot(Sn.stack().values, dist = "norm", plot=pylab)
stats.probplot(Dn.stack().values, dist = "norm", plot=pylab)


# one cell-line S/D, all treatment, all replicates
# groups is S and D


stats.f_oneway(S, D)

S,D,SD = split_surviving_dead(df_log2FC_norm, cell_lines = ["A549", "MCF7", "RKO"])
plt.figure(figsize=(14,6))
stats.probplot(Sn.stack().values, dist = "norm", plot=pylab)
stats.probplot(Dn.stack().values, dist = "norm", plot=pylab)


protein = S.index[0]
S_comp = S[S.index == protein]
D_comp = D[D.index == protein]
anova = stats.f_oneway(S_comp.values[0] ,D_comp.values[0])

len(S)

anova_p = stats.f_oneway(S,D,axis=1)[1]
anova_p = pd.DataFrame(anova_p, columns = ["p"])
anova_p["p"].hist(bins=100)
anova_p = anova_p.sort_values(by = "p")
anova_p["q"] = qvalues(anova_p, pcolname = "p")
anova_p["q"].hist(bins=1000)
anova_p.dropna().q.hist(bins=200)

protein

# two treatments
S_t = S.filter(regex="1_.*")
D_t = D.filter(regex="1_.*")

def kruskal_test(S,D):
    protein_array = []
    kruskal = []
    for protein in list(set(S.index) & set(D.index)):        
        S_val = S[S.index == protein].values[0]
        if np.isnan(S_val).all():
            continue
        D_val = D[D.index == protein].values[0]
        if np.isnan(D_val).all():
            continue
        kruskal_p = stats.kruskal(S_val,D_val,nan_policy = "omit")[1]
        kruskal.append(kruskal_p)
        protein_array.append(protein)
    
    kruskal = pd.DataFrame(kruskal, index = protein_array, columns = ["p"])
    kruskal.sort_values(by = "p", inplace=True)
    kruskal["q"] = qvalues(kruskal, pcolname = "p")
    return kruskal
    


start = time.time()
kruskal_test(S,D)
end = time.time()
print(end-start)
kruskal[kruskal.q > 0.05]
kruskal[kruskal.p > 0.05]



protein

test = stats.kruskal(
S.filter(regex="1_A549*"),
S.filter(regex="2_A549*"),
S.filter(regex="3_A549*"),
S.filter(regex="4_A549*"),
S.filter(regex="5_A549*"),
S.filter(regex="6_A549*"),
S.filter(regex="7_A549*"),
S.filter(regex="8_A549*"),
S.filter(regex="9_A549*"), nan_policy = "omit")


cell_line = "A549"

def kruskal_treatment_groups(S,D, cell_line = "all"):
    s_1 = S.filter(regex=f"1_{cell_line}*")
    s_2 = S.filter(regex=f"2_{cell_line}*")
    s_3 = S.filter(regex=f"3_{cell_line}*")
    s_4 = S.filter(regex=f"4_{cell_line}*")
    s_5 = S.filter(regex=f"5_{cell_line}*")
    s_6 = S.filter(regex=f"6_{cell_line}*")
    s_7 = S.filter(regex=f"7_{cell_line}*")
    s_8 = S.filter(regex=f"8_{cell_line}*")
    s_9 = S.filter(regex=f"9_{cell_line}*")
    d_1 = D.filter(regex=f"1_{cell_line}*")
    d_2 = D.filter(regex=f"2_{cell_line}*")
    d_3 = D.filter(regex=f"3_{cell_line}*")
    d_4 = D.filter(regex=f"4_{cell_line}*")
    d_5 = D.filter(regex=f"5_{cell_line}*")
    d_6 = D.filter(regex=f"6_{cell_line}*")
    d_7 = D.filter(regex=f"7_{cell_line}*")
    d_8 = D.filter(regex=f"8_{cell_line}*")
    d_9 = D.filter(regex=f"9_{cell_line}*")
    
    if cell_line == "all":
        s_1 = S.filter(regex=f"1_.*")
        s_2 = S.filter(regex=f"2_.*")
        s_3 = S.filter(regex=f"3_.*")
        s_4 = S.filter(regex=f"4_.*")
        s_5 = S.filter(regex=f"5_.*")
        s_6 = S.filter(regex=f"6_.*")
        s_7 = S.filter(regex=f"7_.*")
        s_8 = S.filter(regex=f"8_.*")
        s_9 = S.filter(regex=f"9_.*")
        d_1 = D.filter(regex=f"1_.*")
        d_2 = D.filter(regex=f"2_.*")
        d_3 = D.filter(regex=f"3_.*")
        d_4 = D.filter(regex=f"4_.*")
        d_5 = D.filter(regex=f"5_.*")
        d_6 = D.filter(regex=f"6_.*")
        d_7 = D.filter(regex=f"7_.*")
        d_8 = D.filter(regex=f"8_.*")
        d_9 = D.filter(regex=f"9_.*")
    
    protein_array = []
    kruskal_p = []
    for protein in d_1.index:
        try:
            kruskal = stats.kruskal(s_1[s_1.index == protein],
                                      s_2[s_2.index == protein],
                                      s_3[s_3.index == protein],
                                      s_4[s_4.index == protein],
                                      s_5[s_5.index == protein],
                                      s_6[s_6.index == protein],
                                      s_7[s_7.index == protein],
                                      s_8[s_8.index == protein],
                                      s_9[s_9.index == protein],  
                                      d_1[d_1.index == protein],
                                      d_2[d_2.index == protein],
                                      d_3[d_3.index == protein],
                                      d_4[d_4.index == protein],
                                      d_5[d_5.index == protein],
                                      d_6[d_6.index == protein],
                                      d_7[d_7.index == protein],
                                      d_8[d_8.index == protein],
                                      d_9[d_9.index == protein], nan_policy = "omit")
            protein_array.append(protein)
            kruskal_p.append(kruskal[1])
        except:
            continue    
    
    kruskal_df = pd.DataFrame(kruskal_p, index = protein_array, columns = ["p"])
    kruskal_df.sort_values(by="p",inplace=True)
    kruskal_df["q"] = qvalues(kruskal_df, pcolname="p")
    return kruskal_df



kruskal_df.p.hist(bins=100)
kruskal_df.q.hist(bins=100)

s_1[s_1.index == protein]
S.filter(regex="2_A549*"), nan_policy="omit")


.dropna()
D.dropna()


a = np.array([[9.87, 9.03, 6.81],
              [7.18, 8.35, 7.00],
              [8.39, 7.58, 7.68],
              [7.45, 6.33, 9.35],
              [6.41, 7.10, 9.33],
              [8.00, 8.24, 8.44]])
b = np.array([[6.35, 7.30, 7.16],
              [6.65, 6.68, 7.63],
              [5.72, 7.73, 6.72],
              [7.01, 9.19, 7.41],
              [7.75, 7.87, 8.30],
              [6.90, 7.97, 6.97]])
c = np.array([[3.31, 8.77, 1.01],
              [8.25, 3.24, 3.62],
              [6.32, 8.81, 5.19],
              [7.48, 8.83, 8.91],
              [8.59, 6.01, 6.07],
              [3.07, 9.72, 7.48]])
    
    
stats.f_oneway(a,b,c, axis = 1)








anova = pd.DataFrame(stats.f_oneway(Sn, Dn, axis = 1)[0], index = Sn.index, columns = ["p"])
anova.sort_values(by = "p", inplace=True)
anova["q"] = qvalues(anova, pcolname = "p")

anova.p.hist(bins=100)
anova.q.hist(bins=100)

anova = pd.DataFrame(stats.f_oneway(Sn.filter(regex=f"1_A549"), Dn.filter(regex=f"1_A549"), axis = 1)[0], index = Sn.index, columns = ["p"])
anova.sort_values(by = "p", inplace=True)
anova["q"] = qvalues(anova, pcolname = "p")

anova.q.hist(bins=100)


anova[anova.q<0.0]







