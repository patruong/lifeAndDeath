#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:02:25 2021

@author: ptruong
"""

import pandas as pd
import numpy as np



df = pd.read_csv("HBY20Mix.peptides.csv", index_col = 0)
df = pd.read_csv("iPRG.novo.pep.csv", index_col = 0)


np.log(df.iloc[:,0:9]).hist(bins=100)

import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=1, ncols=1)

for ind, col in enumerate(df.columns[1:10]):
    np.log(df[col]).hist(bins=500, alpha = 0.5, histtype="step")


