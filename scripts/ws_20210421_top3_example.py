#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:38:29 2021

@author: ptruong
"""

import pandas as pd
import time

df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")



from scipy.stats import ttest_ind

start = time.time()
for protein, peptides in df_raw.groupby("Leading razor protein"):

  reduced = peptides.copy().filter(regex=("^Reporter intensity corrected [01] A549_D.*"))
  reduced['sum'] = reduced.sum(axis=1)
  case_col = [col for col in reduced if col.startswith("Reporter intensity corrected 0")]
  ctrl_col = [col for col in reduced if col.startswith("Reporter intensity corrected 1")]
  reduced.sort_values(by="sum", ascending=False, inplace=True)
  reduced = reduced.iloc[0:3,:][case_col+ctrl_col]
  testable = reduced.sum(axis=0)
  p = ttest_ind(testable[case_col], testable[ctrl_col])[1]
  if p<0.005:
    print(protein,p)

end = time.time()
print(end-start)



