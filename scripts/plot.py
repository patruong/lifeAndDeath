#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:53:18 2021

@author: ptruong
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer


def get_variables():
    cell_lines = ['A549', 'MCF7', 'RKO']
    states = ['D', 'S']
    replicates = ['Rep1', 'Rep2', 'Rep3']
    return cell_lines, states, replicates


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
    cell_lines, states, replicates = get_variables()
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
    cell_lines, states, replicates = get_variables()
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
    cell_lines, states, replicates = get_variables()
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
    colour = "treatment"
    
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
    
    
    
 

################
# VOLCANO PLOT #
################



def get_significant_proteins(df_log2FC, df_pVals, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05):
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

def volcano_plot(df_log2FC, df_pVals, cell_line, state, treatment, fc_treshold = 1.0, pVal_treshold = 0.05):
    df_volcano = get_significant_proteins(df_log2FC, df_pVals, cell_line, state, treatment, fc_treshold = fc_treshold, pVal_treshold = pVal_treshold)
    b = sns.scatterplot(x="log2FC", y="-log10(p-value)", hue = "sign", data = df_volcano)
    for i in range(np.shape(df_volcano[df_volcano.sign != "non-significant"])[0]):
        vals = df_volcano[df_volcano.sign != "non-significant"].iloc[i,:]
        plt.text(vals["log2FC"], -np.log10(vals["p-value"]), vals.name , 
                 horizontalalignment='left', size='small', color='black', weight='light')
    b.set_xlabel("log2FC",fontsize=4)
    plt.show()

   
    