#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 00:10:06 2021

@author: ptruong
"""


# Mapping function to map column name to seperate columns with cell line, state, treatment or replicate.

col_to_treatment_mapper = lambda x : x.split(" ")[3]
treatment_nomenclature_map_dict = {"0": "control", "1": "8-azaguanine", "2": "Raltitrexed", "3": "Topotecan", "4":"Floxuridine", "5":"Nutlin", "6":"Dasatanib", "7":"Gefitinib", "8":"Vincristine", "9":"Bortezomib"}
col_to_cell_line_mapper = lambda x : x.split(" ")[4].split("_")[0]
col_to_state_mapper = lambda x : x.split(" ")[4].split("_")[1] 
col_to_rep_mapper = lambda x : x.split(" ")[4].split("_")[2]
