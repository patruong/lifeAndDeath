#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:21:23 2021

@author: ptruong
"""




import os
import time

import pandas as pd
import numpy as np

from load import load_peptide_data
from top3 import top3
from get_variables import get_cell_lines_states_replicates

df = load_peptide_data("peptides tryptic.csv", max_PEP = 0.1, max_missed_cleavages = 3)

df = df.iloc[df.index.get_level_values("Unique (Groups)") == 'yes', :] # Take only peptides with unique protein groups

cell_lines,states, replicates = get_cell_lines_states_replicates()
treatments = [str(i) for i in range(10)]

# Take Top3 of these 

df.iloc[:, df.columns.get_level_values("sample") == "A549_D_1"]




