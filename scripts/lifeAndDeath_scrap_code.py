#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 19:13:44 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

data = pd.read_csv("proteinGroups_tryptic_reporterIntensities_noNorm.csv", index_col = 0)
meta = pd.read_csv("proteinGroups_tryptic_sampleMetadata_noNorm.csv", index_col = 0)

data_col = data.columns
tuples = list(zip(data_col.values, meta.treatment.values, meta.cell_line.values, meta.state.values, meta.replicate.values))
index = pd.MultiIndex.from_tuples(tuples, names=['data', 'treatment', "cell_line", "state", "repl"])


df = pd.DataFrame(data.values.T, index = index, columns = data.index)


df_dead = df.iloc[df.index.get_level_values("state") == "D"]
df_alive = df.iloc[df.index.get_level_values("state") == "S"]


# Assumption replace -inf with 0
array_logFc = np.log2(df_alive).replace(-np.inf, 0).values - np.log2(df_dead).replace(-np.inf, 0).values
array_logFc_ = np.log2(df_dead).replace(-np.inf, 0).values - np.log2(df_alive).replace(-np.inf, 0).values


array_logFc_ = np.log2(df_dead).replace(-np.inf, 0).values - np.log2(df_alive).replace(-np.inf, 0).values


idx = df_dead.index
col = df_dead.columns
df_logFc = pd.DataFrame(array_logFc, index = idx, columns = col) # log2Fc alive - dead

# columnwise differential expression


# Rethink how to compute the logFC -> should be from base-line
# T-test should be between life and dead cells...


df_control = df.iloc[df.index.get_level_values("treatment") == "Control"]
df_treatment1 = df.iloc[df.index.get_level_values("treatment") == "Control"]
df_treatment2 = df.iloc[df.index.get_level_values("treatment") == "Nutlin"]
df_treatment1 = df.iloc[df.index.get_level_values("treatment") == "Topotecan"]

"""
array_logFc = np.log2(df_control).replace(-np.inf, 0).values - np.log2(df_treatment1).replace(-np.inf, 0).values
idx = df_control.index
col = df_control.columns
df_logFc = pd.DataFrame(array_logFc, index = idx, columns = col) # log2Fc control vs treatment


array_logFc = np.log2(df_control).replace(-np.inf, 0).values - np.log2(df_treatment2).replace(-np.inf, 0).values
idx = df_control.index
col = df_control.columns
df_logFc2 = pd.DataFrame(array_logFc, index = idx, columns = col) # log2Fc control vs treatment

df_logFc.append(df_logFc2)
"""

df_control = df.iloc[df.index.get_level_values("treatment") == "Control"]
df_logFc = pd.DataFrame()
for treatment in meta.treatment.unique():
    df_treatment = df.iloc[df.index.get_level_values("treatment") == treatment]
    #Note imputing -np.inf values with 0
    array_logFc = np.log2(df_control).replace(-np.inf, 0).values - np.log2(df_treatment).replace(-np.inf, 0).values
    idx = df_control.index
    col = df_control.columns
    df_tmp = pd.DataFrame(array_logFc, index = idx, columns = col) # log2Fc control vs treatment
    df_logFc = df_logFc.append(df_tmp)
    
    

# Ask P about travel to amsterdam! 

c = pd.DataFrame()
a = pd.DataFrame([[1,2],[2,3]], columns = ["A", "B"])
b = pd.DataFrame([[3,2],[2,3]], columns = ["A", "B"])

a.append(b)

meta.treatment.unique()

#pd.DataFrame(array_logFc).isnull().sum().sum()
#pd.DataFrame(array_logFc_).isnull().sum().sum()


pd.DataFrame(array_logFc, index = index, columns = data.index)







