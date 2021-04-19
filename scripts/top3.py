#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 18:15:12 2021

@author: ptruong
"""

import time
import numpy as np
import pandas as pd
import scipy 

#########################
# Protein summarization #
#########################
        
#df = data_aitchison

#df = df_int

def top3(df, output_name = "top3.protein.tsv"):
    print("Running top3 summarization")
    print("NOTE: runtime for this on planck is about 25000seconds or 7 hours.")
    df = pd.DataFrame(df.values, index = df.index, columns = df.columns.get_level_values("experiment"))
    start = time.time()
    protein_quants = []
    experiment_array = []
    for col_i in range(len(df.columns)):
        protein_quant_array = []
        for protein in df.index.unique():
            top3_peptides = df[df.index == protein].iloc[:,col_i].nlargest(3)
            if len(top3_peptides)>=2:
                proteinQuant = top3_peptides.mean()
            else:
                proteinQuant = np.nan
            protein_quant_array.append(proteinQuant)
        print(time.time()-start)
        protein_quants.append(protein_quant_array)
        experiment_array.append(top3_peptides.name)
    
    df_protQuant = pd.DataFrame(protein_quants, index = experiment_array, columns = df.index.unique()).T
    df_protQuant.to_csv(output_name, sep = "\t")
    end = time.time()
    print(end-start)
    print("top3 finished")
    return df_protQuant

###########################
# Protein sum col to mIdx #
###########################

def protSum_col_to_mIdx(df_int):
    treatment_mapper = lambda x : x.split("_")[0]
    specie_mapper = lambda x : x.split("_")[1]
    state_mapper = lambda x : x.split("_")[2]
    replicate_mapper = lambda x : x.split("_")[3]
    batch_mapper = lambda x : x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[3]
    experiment_mapper = lambda x : x.split("_")[0] + "_" + x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[3]
    sample_mapper = lambda x : x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[0]

    def map_values_to_value_list(value_list, values):
        """
        e.g. test below lists to see results.
        value_list = ["a", "b", "c"]
        values = ["c","c","a"]
        """
        return [value_list.index(x) for x in values]
    
    specie_col = df_int.columns.map(specie_mapper)
    specie_list = list(specie_col.unique())
    specie_code = map_values_to_value_list(specie_list, specie_col)
    
    state_col = df_int.columns.map(state_mapper)
    state_list = list(state_col.unique())
    state_code = map_values_to_value_list(state_list, state_col)
    
    replicate_col = df_int.columns.map(replicate_mapper)
    replicate_list = list(replicate_col.unique())
    replicate_code = map_values_to_value_list(replicate_list, replicate_col)
    
    treatment_col = df_int.columns.map(treatment_mapper)
    treatment_list = list(treatment_col.unique())
    treatment_code = map_values_to_value_list(treatment_list, treatment_col)
    
    batch_col = df_int.columns.map(batch_mapper)
    batch_list = list(batch_col.unique())
    batch_code = map_values_to_value_list(batch_list, batch_col)
    
    sample_col = df_int.columns.map(sample_mapper)
    sample_list = list(sample_col.unique())
    sample_code = map_values_to_value_list(sample_list, sample_col)

    experiment_col = df_int.columns.map(experiment_mapper)
    experiment_list = list(experiment_col.unique())
    experiment_code = map_values_to_value_list(experiment_list, experiment_col)

    
    midx = pd.MultiIndex(levels=[specie_list, state_list, treatment_list, replicate_list, batch_list, sample_list, experiment_list], 
                         codes=[specie_code, state_code, treatment_code, replicate_code, batch_code, sample_code, experiment_code],
                         names=["cell_line", "state", "treatment", "replicate", "batch", "sample", "experiment"])
    return midx


def protSum_intensities_to_midx_df(df_int):
    midx = protSum_col_to_mIdx(df_int)
    df_int.columns = midx
    return df_int



def aggregate_protein_quantity(df_protQuant):
    """
    Get the mean aggregate of protein quantites.
    """
    df_quant = pd.DataFrame()
    for sample in df_protQuant.columns.get_level_values("sample").unique():
        df_sample = df_protQuant.iloc[:,df_protQuant.columns.get_level_values("sample") == sample]
        #count = df_sample.count(axis=1) > 1 #sample treshold
        #df_sample = (df_sample.T*count.values).T
        df_sample_quant = pd.DataFrame(df_sample.mean(axis=1), columns = [sample])
        df_quant = pd.concat([df_quant, df_sample_quant], axis = 1)
    return df_quant.replace(0, np.nan)

def get_p_matrix(df_protQuant):
    df_pVals = pd.DataFrame()
    df_tStats = pd.DataFrame()
    for treatment in df_protQuant.columns.get_level_values("treatment").unique()[1:]:
        df_control = df_protQuant.iloc[:,df_protQuant.columns.get_level_values("treatment") == "0"]    
        for cell_line in df_protQuant.columns.get_level_values("cell_line").unique():
            for state in df_protQuant.columns.get_level_values("state").unique():
                sample = cell_line + "_" + state + "_" + treatment
                df_sample = df_protQuant.iloc[:, df_protQuant.columns.get_level_values("sample") == sample]
                df_sample_control = df_control[cell_line][state]
                tStats = pd.DataFrame(scipy.stats.ttest_ind(df_sample,df_sample_control, axis = 1)[0], columns = [sample], index = df_sample.index)
                pVals = pd.DataFrame(scipy.stats.ttest_ind(df_sample,df_sample_control, axis = 1)[1], columns = [sample], index = df_sample.index)
                df_tStats = pd.concat([df_tStats, tStats], axis = 1)
                df_pVals = pd.concat([df_pVals, pVals], axis = 1)
    return df_pVals, df_tStats




