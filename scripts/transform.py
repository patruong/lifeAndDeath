#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:51:54 2021

@author: ptruong
"""
import time

import numpy as np
import pandas as pd

from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def get_variables():
    cell_lines = ['A549', 'MCF7', 'RKO']
    states = ['D', 'S']
    replicates = ['Rep1', 'Rep2', 'Rep3']
    return cell_lines, states, replicates

def aitchison_transform_part(df, use_multiplicative_replacement = True):
    """
    Aitchison tranformation on df with all columns belonging to same batch.
    
    df should consist of all samples tagged together in one channel (i.e. A549_S_rep1 etc.)
    """
    if use_multiplicative_replacement == True:
        df_aitchison = multiplicative_replacement(df)
    else:
        df_aitchison = closure(df)
    df_idx = df.index
    df_col = df.columns
    df_aitchison = pd.DataFrame(df_aitchison, index = df_idx, columns = df_col)
    return df_aitchison

def aitchison_transform(df_int, use_multiplicative_replacement = True):
    cell_lines, states, replicates = get_variables()
    df_aitchison = pd.DataFrame()
    
    start = time.time()
    for cell_line in cell_lines:
        for state in states:
            for rep in replicates:
                df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == ("_".join([cell_line, state, rep]))]
                df_part = df_part[(df_part.T != 0).any()]
                df_part = aitchison_transform_part(df_part, use_multiplicative_replacement)
                df_aitchison = pd.concat([df_aitchison, df_part], axis = 1 )
                print(time.time()-start)
    end=time.time()
    print(end-start)
    return df_aitchison

   
def norm_SL(df_int):
    target = df_int.sum().mean()
    norm_fac = target/df_int.sum()
    df_norm = df_int*norm_fac
    return df_norm


def calcNormFactors(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
      r_from_pd_df = ro.conversion.py2rpy(df)
    
    edgeR = importr("edgeR")
    raw_tmm = edgeR.calcNormFactors(r_from_pd_df)
    
    with localconverter(ro.default_converter + pandas2ri.converter):
      pd_from_r_df = ro.conversion.rpy2py(raw_tmm)
      
    raw_tmm_pd = pd.DataFrame(pd_from_r_df, index =df.columns).T.values
    
    sl_tmm = raw_tmm_pd
    sl_tmm = pd.Series(sl_tmm[0], df.columns)

    return sl_tmm



def irs_norm(df_norm):
    irs = pd.DataFrame()
    for batch in df_norm.columns.get_level_values("batch").unique():
        irs_sum = df_norm.iloc[:,df_norm.columns.get_level_values("batch") == batch].sum(axis = 1)
        irs_sum = pd.DataFrame(irs_sum.values, columns = ["sum_" + batch], index = df_norm.index)
        irs = pd.concat([irs, irs_sum], axis = 1)
    
    irs_geoAvg = np.exp(np.log(df_norm.replace(0, np.nan)).mean(axis=1))
    irs_geoAvg = pd.DataFrame(irs_geoAvg, columns = ["geoAvg"])
    irs = pd.concat([irs, irs_geoAvg], axis = 1)
    irs = irs.replace(0,np.nan)
    
    i = 0
    for batch in df_norm.columns.get_level_values("batch").unique():
        irs_fac = irs.geoAvg/irs.iloc[:,i]
        irs_fac = pd.DataFrame(irs_fac.values, columns = ["fac_"+batch], index = df_norm.index)
        irs = pd.concat([irs, irs_fac], axis = 1)    
        i+=1
    
    df_irs = pd.DataFrame()
    for batch in df_norm.columns.get_level_values("batch").unique():
        exp_sl = df_norm.iloc[:,df_norm.columns.get_level_values("batch") == batch]
        irs_fac = pd.DataFrame(irs["fac_" + batch])
        
        data_irs = irs_fac.values*exp_sl
        df_irs = pd.concat([df_irs, data_irs], axis = 1)
    return df_irs





