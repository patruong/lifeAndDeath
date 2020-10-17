#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 21:02:29 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

from uniprot import *

df = pd.read_csv("melted_treshold.csv", sep ="\t")

df.Leading_razor_protein[0]
#df.Proteins[0]
df_sub = df[0:200]

import time

start_time = time.time()
tmp_time = start_time
go_list = []
for i in df_sub.Leading_razor_protein:
    go_list.append(go_class(i))
    proc_time = time.time()
    print("From start: " + str(proc_time - start_time))
    print("Iter time: " + str(proc_time - tmp_time))
    tmp_time = proc_time
    


start_time=time.time()
df_go = df_sub.Leading_razor_protein.apply(go_class)
end_time=time.time()
print(end_time - start_time)







