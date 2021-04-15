#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:29:20 2021

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
from get_variables import get_cell_line_states_replicates_from_reporter_intensity_cols
from midx import col_to_mIdx, intensities_to_midx_df, diffacto_col_to_mIdx, diffacto_to_midx_df
from pd_functions import drop_zero_row
from transform import aitchison_transform_part, aitchison_transform, norm_SL, calcNormFactors, irs_norm
from plot import plot_kde, kde_matrix_plot_all_channels, plot_kde_batch, kde_matrix_plot_batch, plot_intensity_boxplot, plot_diffacto_pca, plot_diffacto_pca_cell_line, pca_plot_ax, get_significant_proteins, volcano_plot 
from top3 import top3, protSum_col_to_mIdx, protSum_intensities_to_midx_df, aggregate_protein_quantity, get_p_matrix
from de_analysis import get_log2FC_regulation, get_log2FC

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

df_int = df[reporter_intensity_corrected_cols]
df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
df_int = drop_zero_row(df_int) # dropped last zero-row
midx = col_to_mIdx(df_int)

#######
# Raw #
#######
df_int = intensities_to_midx_df(df_int)
df_signal = (df_int != 0) # create a boolean matrix for true and false.

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


##############
# AITCHISOLN #
##############

data_aitchison = aitchison_transform(df_irs_tmm, use_multiplicative_replacement = True)
kde_matrix_plot_batch(data_aitchison, log = False, suptitle = "Aitchison(IRS-TMM Peptide) (Channel intensities together)", legend_size = 6)

data_raw_aitchison = aitchison_transform(df_int.replace(0,np.nan), use_multiplicative_replacement = True)
kde_matrix_plot_batch(data_raw_aitchison, log = False, suptitle = "Aitchison(Raw Peptide) (Channel intensities together)", legend_size = 6)


###########################
# AITCHISON PEPTIDE PCA  ##
###########################

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


#########
# TOP3 ##
#########

df_protQuant = top3(data_raw_aitchison, output_name = "PROT_QUANT_DF_aitchison_no_norm.csv")
df_protQuant = pd.read_csv("PROT_QUANT_DF.csv", sep="\t", index_col=0, header=[0,1,2,3,4,5,6])
df_protQuant = protSum_intensities_to_midx_df(df_protQuant)


kde_matrix_plot_batch(df_protQuant, log=False, suptitle = "Raw->SL->IRS->TMM-Aitchison->Top3") #looks ok?

#############
# TOP3 PCA ##
#############

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

##################
# Aggregate TOP3 #
###################

df_quant= aggregate_protein_quantity(df_protQuant)
df_quant = diffacto_to_midx_df(df_quant, protein_index=False) 
df_pVals, df_tStats = get_p_matrix(df_protQuant)


plt.subplot(2, 2, 1)
plot_diffacto_pca(df_quant) # JUMP DOWN TO DIFFACTO SECTION TO GET PLOTTING FUNCTION....
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (all)")
plt.subplot(2, 2, 2)
plot_diffacto_pca_cell_line(df_quant.A549) # Good !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (A549)")
plt.subplot(2, 2, 3)
plot_diffacto_pca_cell_line(df_quant.RKO) # Good  !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (RKO)")
plt.subplot(2, 2, 4)
plot_diffacto_pca_cell_line(df_quant.MCF7) # Good !!!!!!
plt.title("Raw->SL-IRS->TMM->Aitchison->Top3 PCA (MCF7)")

###########
# DE TOP3 #
###########


df_log2FC, df_pVals, df_tStats = get_log2FC(df_protQuant)


cell_line = "A549"
state = "D"
treatment = "1"
df_significant = get_significant_proteins(df_log2FC, df_pVals, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05)


cell_line = "A549"
state = "S"
fc_tresh = 0.7
pVal_tresh = 0.05

for treatment in range(1,10):
    i = treatment
    plt.subplot(3,3,i)
    df_significant = get_significant_proteins(df_log2FC, df_pVals, cell_line, state, str(treatment), fc_treshold = fc_tresh, pVal_treshold = pVal_tresh)
    n_up_reg = (df_significant.sign == "up-regulated").sum()
    n_down_reg = (df_significant.sign == "down-regulated").sum()
    volcano_plot(df_log2FC, df_pVals, cell_line, state, str(treatment), fc_treshold = fc_tresh, pVal_treshold = pVal_tresh)
    plt.title("treatment: " + str(treatment) +  " up: " + str(n_up_reg) + " down: " + str(n_down_reg))    
#plt.tight_layout()
plt.suptitle(cell_line + "_" + state + " "+ "FC: " + str(fc_tresh) + " p-value: " + str(pVal_tresh))




############
# Diffacto #
############

# check we_20210412 for data_aitchison to diffact input conversion

os.chdir("/home/ptruong/git/lifeAndDeath/scripts/diffacto")
diffacto = pd.read_csv("diffacto_output.protein.txt", sep = "\t")
diffacto_FDR = pd.read_csv("diffacto_output.protein.FDR", sep = "\t")
diffacto["MCFDR"] = diffacto_FDR["MCFDR"]
#diffacto = diffacto[diffacto["MCFDR"] < 0.05] #treshold on monte-carlo FDR
#diffacto = diffacto[diffacto["S/N"] > 0] # treshold on signal-to-noise ratio

df_diffacto = diffacto_to_midx_df(diffacto)



plot_diffacto_pca(df_diffacto)
plt.title("PCA (Aitchison -> Diffacto)")

plot_diffacto_pca_cell_line(df_diffacto.A549)
plt.title("PCA (Aitchison -> Diffacto)")