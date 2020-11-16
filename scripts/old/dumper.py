#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 23:41:28 2020

@author: ptruong
"""


from __future__ import print_function
import time
import numpy as np
import pandas as pd

from preprocessor import *
from stats import *

def get_sample_metadata(data):
    """
    Get sample metadata from create_full_data() function.
    """
    data_values, target_drugs, cell_lines, states, replicates = split_data(data)
    treatment_list = ["Control",
                    "8-zaguanine",
                    "Raltitrexed",
                    "Topotecan",
                    "Floxuridine",
                    "Nutlin",
                    "Dasatinib",
                    "Gefitinib",
                    "Vincristine",
                    "Bortezomib"]
    treatment_number = list(range(0,10))
    treatment_dict = dict(zip(treatment_number, treatment_list))
    mapper = lambda t: treatment_dict[t]
    vfunc = np.vectorize(mapper)
    treatments = vfunc(target_drugs)
    df = pd.DataFrame(np.array([treatments, cell_lines, states, replicates]).T, index = data.columns, columns=["treatment", "cell_line", "state", "replicate"])
    return df


def rename_cols(data):
    """
    Renames columns of the data
    """
    treatment_list = ["Control",
                    "8-zaguanine",
                    "Raltitrexed",
                    "Topotecan",
                    "Floxuridine",
                    "Nutlin",
                    "Dasatinib",
                    "Gefitinib",
                    "Vincristine",
                    "Bortezomib"]
    treatment_number = list(range(0,10))
    treatment_dict = dict(zip(treatment_number, treatment_list))
    col_names = []
    for i in data.columns:
        col = i[29:-5]
        col_name = treatment_dict[int(col.split()[0])] + "_" + col.split()[1]
        col_names.append(col_name)
    renamed_cols = dict(zip(data.columns, col_names))
    data = data.rename(columns = renamed_cols)
    return data


if __name__ == "__main__":    
    filename = "proteinGroups tryptic.txt"
    data = create_full_data(data_file = filename, treshold = 1)
#    data = data.fillna(0) # Input missing values to zero
    data = normalize(data)
    data.fillna(0)
    sample_meta_data = get_sample_metadata(data)
    data_new_cols = rename_cols(data) #just renames cols
    
    sample_meta_data.to_csv("proteinGroups_tryptic_sampleMetadata_noNorm.csv", sep = ",")
    data_new_cols.to_csv("proteinGroups_tryptic_reporterIntensities_noNorm.csv", sep = ",")

