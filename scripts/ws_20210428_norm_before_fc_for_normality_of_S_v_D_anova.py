#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 17:41:47 2021

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
from de_analysis import get_log2FC_regulation, get_log2FC
from q_value import qvalues

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



#######################
# DrOP BAD replicates #
#######################


kde_matrix_plot_all_channels(intensities_to_midx_df(df_int), log = True, suptitle ="log2(raw)")
kde_matrix_plot_batch(intensities_to_midx_df(df_int), log = True, suptitle = "log2(RawPeptide) (Channel intensities together)", legend_size = 6)


# We drop columns that appear off in one replicate but not other... for example t9 MCF7 is higher in all samples... this is ok!
bad_rep = ["Reporter intensity corrected 0 A549_D_Rep1", "Reporter intensity corrected 3 RKO_D_Rep1", "Reporter intensity corrected 3 RKO_D_Rep3"]
df_int = df_int.drop(bad_rep, axis = 1)

midx = col_to_mIdx(df_int)

#######
# Raw #
#######
df_int = intensities_to_midx_df(df_int)
df_signal = (df_int != 0) # create a boolean matrix for true and false.

kde_matrix_plot_all_channels(df_int, log = True, suptitle ="log2(norm)")
kde_matrix_plot_batch(df_int, log = True, suptitle = "log2(RawPeptide) (Channel intensities together)", legend_size = 6)

########
# NORM #
########
    
df_norm = norm_SL(df_int)
kde_matrix_plot_all_channels(df_norm, log = False, suptitle ="log2(norm)")
kde_matrix_plot_batch(df_norm, log = False, suptitle = "log2(Norm intensity) (Channel intensities together)", legend_size = 6)


#######
# TMM #
#######

sl_tmm = calcNormFactors(df_norm)
df_tmm = df_norm/(sl_tmm)

kde_matrix_plot_all_channels(df_tmm, log = True, suptitle ="log2(TMM)")
kde_matrix_plot_batch(df_tmm, log = True, suptitle = "log2(TMM) (Channel intensities together)", legend_size = 6)


#######
# IRS #
#######

df_irs = irs_norm(df_norm)

kde_matrix_plot_all_channels(df_irs, log = True, suptitle ="log2(IRS)")
kde_matrix_plot_batch(df_irs, log = True, suptitle = "log2(IRS) (Channel intensities together)", legend_size = 6)


#############
# IRS - TMM #
#############

irs_tmm = calcNormFactors(df_irs.fillna(0))
df_irs_tmm = df_irs/(irs_tmm)

kde_matrix_plot_all_channels(df_irs_tmm, log = True, suptitle ="log2(IRS-TMM)")
kde_matrix_plot_batch(df_irs_tmm, log = True, suptitle = "log2(IRS-TMM peptide intensities) (Channel intensities together)", legend_size = 6)

###########################
# FC control vs treatment #
###########################

cell_lines = [ "A549", "MCF7", "RKO"]
states = ["S", "D"]
replicates = ["Rep1", "Rep2", "Rep3"]
treatments = [i for i in range(10)]

df_int = df_int.replace(0, np.nan)
df_log2 = np.log2(df_int)

cell_line = "A549"
state = "S"
rep = "Rep1"

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

S = df_fc["MCF7"].iloc[:, (df_fc["MCF7"].columns.get_level_values("state") == "S")]
D = df_fc["MCF7"].iloc[:, df_fc["MCF7"].columns.get_level_values("state") == "D"]
S = pd.DataFrame(S.values, index = S.index, columns = S.columns.get_level_values("experiment"))
D = pd.DataFrame(D.values, index = D.index, columns = D.columns.get_level_values("experiment"))
col_mapper = lambda x: x.split("_")[0] + "_" + x.split("_")[1] + "_" + x.split("_")[-1] 
S = S.rename(columns = col_mapper)
D = D.rename(columns = col_mapper)
SD = pd.concat([S.stack() ,D.stack()], axis = 1).rename(columns={0:"S", 1:"D"})

sns.scatterplot(data = SD, x = "S", y = "D")
plt.title("FC S vs D ")


#######################
# ASSUMPTION CHECKING #
#######################

# same thing but for summed protein, perform t-test to filter away trash
def shapiro(df):
    shapiro = pd.DataFrame(df.apply(stats.shapiro, axis = 1).tolist(), index = df.index)
    shapiro["q"] = qvalues(shapiro, pcolname="pvalue") 
    return shapiro

shapiro_S = shapiro(S) #Shapiro p > 0.05 indicates normality
shapiro_D = shapiro(D)
shapiro_S.pvalue.hist()
shapiro_D.pvalue.hist()
S_idx = shapiro_S[shapiro_S.q > 0.05].index
D_idx = shapiro_D[shapiro_D.q > 0.05].index
idx = list(set(S_idx) & set(D_idx))


stats.probplot(S.stack().values, dist = "norm", plot=pylab)
stats.probplot(D.stack().values, dist = "norm", plot=pylab)

measurements = np.random.normal(loc = 20, scale = 5, size=100)   
stats.probplot(measurements, dist="norm", plot=pylab)
pylab.show()





