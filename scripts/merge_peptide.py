#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:27:58 2020

@author: ptruong
"""

import os
import time

import pandas as pd
import numpy as np

def merge_peptide_dfs():
#    file_prefix = "peptide_tryptic_logFC2_i"
    
    files = sorted(os.listdir())
    res = pd.DataFrame()
    start_time = time.time()
    tmp_time = start_time
    #n_files = int(len(files)/2)
    for i in range(180):
        filename = files[i]
        print("Reading in " + filename + " : "+ str(i) + "/" + str(len(files)))
        df = pd.read_csv(filename, sep = "\t")
        df_PEP_treshold = df[df["PEP"] < 0.01] # PEP treshold
        df_PEP_MissedCleaveges_treshold = df_PEP_treshold[df_PEP_treshold["Missed_Cleaveges"] < 2] #missed lceaves treshold
        del df_PEP_treshold
        df = df_PEP_MissedCleaveges_treshold
        del df_PEP_MissedCleaveges_treshold
        df = df.drop(["Missed_Cleaveges", "Charges", "Score", "PEP", "Gene_names"], axis = 1) # Drop unnessasary columns
        res = res.append(df)
        proc_time = time.time()
        print("From start:" + str(proc_time - start_time))
        print("Iter time:" + str(proc_time - tmp_time))
        tmp_time = proc_time
        
    print("Process finished: " + str(time.time() - start_time))
    res.to_csv("peptide_melted.csv", sep = "\t", index = False)
    print("DONE!")
    
    
    
