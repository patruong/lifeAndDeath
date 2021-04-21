#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 23:45:52 2021

@author: ptruong
"""


"""
Abstract:
    
    This script is to analyze the log2FC for unprocessed data.
    The theory is that since we are looking within each channel we should
    not need to normalize the data. 
    
    We are looking for proteins that behave different for different treatment.
    In other wrods, proteins that are up-regulated in one treatment, but
    down for another. It should behave similarily across technical replicates.
    We are in other words looking for difference in difference.
    
    
    We are looking into using MANOVA for this.
"""


import os

import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from top3 import protSum_intensities_to_midx_df, aggregate_protein_quantity, get_p_matrix

os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")

df = pd.read_csv("PROT_QUANT_DF_raw_PEP_0.1_missed_cleaves_3.csv", sep = "\t", index_col = 0)
df = protSum_intensities_to_midx_df(df)


df.columns.get_level_values("experiment")






