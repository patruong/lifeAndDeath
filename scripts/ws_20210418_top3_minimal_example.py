#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 23:26:00 2021

@author: ptruong
"""

import os

import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from load import load_peptide_data
from top3 import top3

os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")

df = load_peptide_data("peptides tryptic.csv", max_PEP = 0.1, max_missed_cleavages = 3)
df_protQuadfnt = top3(df, output_name = "PROT_QUANT_DF_raw_PEP_0.1_missed_cleaves_3.csv")


df_ = pd.read_csv("PROT_QUANT_DF_raw_PEP_0.1_missed_cleaves_3.csv", sep = "\t")


