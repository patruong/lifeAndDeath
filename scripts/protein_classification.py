#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 21:02:29 2020

@author: ptruong
"""

import pandas as pd
import numpy as np
from multiprocessing import  Pool


from uniprot import *

df = pd.read_csv("../data/melted_treshold.csv", sep ="\t")

df.Leading_razor_protein[0]
#df.Proteins[0]
df_sub = df[0:1000]

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



def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df
start_time=time.time()
df_go = parallelize_dataframe(df_sub.Leading_razor_protein, go_class)
end_time=time.time()
print(end_time - start_time)


start_time=time.time()
var = pool.map(go_class, df_sub.Leading_razor_protein)
end_time=time.time()
print(end_time - start_time)
