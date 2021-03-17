# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


from functools import reduce

import pandas as pd 
import numpy as np

from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

df = pd.read_csv("peptides tryptic.csv", sep = "\t")
rep_list = list(pd.read_csv("reporter_intensity_list.csv", sep = "\t").columns)

def get_base_cols():
    cols = ["Sequence", "Missed cleavages", "Proteins", "Leading razor protein",
       "Unique (Groups)", "Unique (Proteins)", "PEP", "Score"] 
    return cols

def get_all_intensities():
    rep_list_corrected = []
    for i in rep_list:
        if i.split(" ")[2] == "corrected":
            rep_list_corrected.append(i)
    return rep_list_corrected

def select_rep_state_intensities(rep, state):
    rep = rep
    state = state
    
    rep_list_corrected = []
    for i in rep_list:
        if i.split(" ")[2] == "corrected":
            if i.split(" ")[4].split("_")[2] == ("Rep" + str(rep)):
                if i.split(" ")[4].split("_")[1] == state:
                    rep_list_corrected.append(i)
    return rep_list_corrected

df_base = df[get_base_cols()]

# 6 Experiments, one for each replicate and state

def drop_zero_rows(df):
    df = df[(df.T != 0).any()]
    return df

def clr_on_subset(df_subset):
    df_subset = drop_zero_rows(df_subset)
    df_subset = multiplicative_replacement(df_subset)
    df_subset = clr(df_subset)
    return df_subset
    
def preprocess_df(df, rep, state):
    """
    Aitchi transformed subset of data.
    """
    df_subset = df[select_rep_state_intensities(rep, state)]
    cols = df_subset.columns
    df_subset = drop_zero_rows(df_subset) #index should be the same as protein/peptides
    index = df_subset.index
    df_subset = multiplicative_replacement(df_subset)
    df_subset = clr(df_subset)
    df_subset = pd.DataFrame(df_subset, index = index, columns=cols)
    return df_subset

### subset without clr.

s1 = df[select_rep_state_intensities(1, "S")]
s2 = df[select_rep_state_intensities(2, "S")]
s3 = df[select_rep_state_intensities(3, "S")]
d1 = df[select_rep_state_intensities(1, "D")]
d2 = df[select_rep_state_intensities(2, "D")]
d3 = df[select_rep_state_intensities(3, "D")]

### log transform 

delta = 0.00000000000001

s1 = np.log2(s1.replace([0], delta))
s2 = np.log2(s2.replace([0], delta))
s3 = np.log2(s3.replace([0], delta))
d1 = np.log2(s1.replace([0], delta))
d1 = d1.fillna(0)
d2 = np.log2(d2.replace([0], delta))
d3 = np.log2(d3.replace([0], delta))

### standardized without clr

from sklearn import preprocessing

scaler = preprocessing.StandardScaler()

def standardize_df(df):
    names = df.columns
    scaled_df = scaler.fit_transform(df)
    scaled_df = pd.DataFrame(scaled_df, columns=names)
    return scaled_df

s1 = standardize_df(s1)
s2 = standardize_df(s2)
s3 = standardize_df(s3)
d1 = standardize_df(d1)
d2 = standardize_df(d2)
d3 = standardize_df(d3)

def normalize_df(df):
    names = df.columns
    normalized_df = preprocessing.normalize(df)
    normalized_df = pd.DataFrame(normalized_df, columns=names)
    return normalized_df    

s1 = normalize_df(s1)
s2 = normalize_df(s2)
s3 = normalize_df(s3)
d1 = normalize_df(d1)
d2 = normalize_df(d2)
d3 = normalize_df(d3)

#### CLR transform - THIS IS AITCHISON DOES NOT NEED ANY OF ABOVE NORMALIZATION OR STANDARDIZATION OR LOG
# IS ALREADY LOG IN CLR AND NORMALIZED.
    
s1 = preprocess_df(df, 1, "S")
s2 = preprocess_df(df, 2, "S")
s3 = preprocess_df(df, 3, "S")
d1 = preprocess_df(df, 1, "D")
d2 = preprocess_df(df, 2, "D")
d3 = preprocess_df(df, 3, "D")

s1.merge

data_frames = [df_base, s1, s2, s3, d1, d2, d3]
df_merged = reduce(lambda  left,right: pd.merge(left,right,left_index=True, right_index=True,
                                            how='outer'), data_frames)

# last step, check so that len(s1) is the same for all subsets after merging

############################
# Check for batch-effect ###
############################
import matplotlib.pyplot as plt
    

def plot_histogram_all():
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.hist(s1.stack(), bins = 500, alpha = 0.5, label = "Surviving repl 1")
    ax.hist(s2.stack(), bins = 500, alpha = 0.5, label = "Surviving repl 2")
    ax.hist(s3.stack(), bins = 500, alpha = 0.5, label = "Surviving repl 3")
    ax.hist(d1.stack(), bins = 500, alpha = 0.5, label = "Dead repl 1")
    ax.hist(d2.stack(), bins = 500, alpha = 0.5, label = "Dead repl 2")
    ax.hist(d3.stack(), bins = 500, alpha = 0.5, label = "Dead repl 3")
    
    #ax.set_xlim([-3, 4])
    plt.legend()
    plt.show()


def plot_histogram_subset(df, label):
    fig, ax = plt.subplots(figsize=(12, 6))    
    ax.hist(df.stack(), bins = 500, alpha = 0.5, label = label)
    #ax.set_xlim([-3, 4])
    plt.legend()
    plt.show()

plot_histogram_all()
plot_histogram_subset(s1, "Surviving repl 1")
plot_histogram_subset(s2, "Surviving repl 2")
plot_histogram_subset(s3, "Surviving repl 3")
plot_histogram_subset(d1, "Dead repl 1")
plot_histogram_subset(d2, "Dead repl 2")
plot_histogram_subset(d3, "Dead repl 3")



df_tresholded = df_merged[df_merged["PEP"] < 0.01]

# THIS IS GOOD and what we work with...

############

import numpy as np
import skbio.stats.composition
from skbio.stats.composition import perturb

otus = np.array([1./3, 1./3., 1./3])
antibiotic = np.array([1./2, 1./2, 1])
perturb(otus, antibiotic)


from skbio.stats.composition import closure
X = np.array([[2, 2, 6], [4, 4, 2]])
closure(X)

from skbio.stats.composition import clr
x = np.array([.1, .3, .4, .2])
clr(x)
clr(X)


