#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 23:08:45 2021

@author: ptruong
"""

import pandas as pd
import numpy as np

from get_columns import get_all_reporter_intensity_correct
from pd_functions import drop_zero_row
from midx import col_to_mIdx, intensities_to_midx_df

def load_peptide_data(filename, max_PEP = 0.01, max_missed_cleavages = 3):
    df_raw = pd.read_csv(filename, sep = "\t")
    
    reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()
    
    
    # filter df
    df = df_raw[df_raw.PEP < max_PEP] # 5% PEP removed 9578 peptides
    df = df[df["Missed cleavages"] < max_missed_cleavages] # 0 removed
    #df = df.set_index("Leading razor protein")
    df = df.set_index(['Proteins', "Leading razor protein", 'Unique (Proteins)', 'Unique (Groups)', 'PEP', "Score"])
    
    df_int = df[reporter_intensity_corrected_cols]
    df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
    df_int = drop_zero_row(df_int) # dropped last zero-row
    
    #######
    # Raw #
    #######
    df_int = intensities_to_midx_df(df_int)
    return df_int







