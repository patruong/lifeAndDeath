#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:32:35 2020

@author: ptruong
"""


import pandas as pd
import numpy as np



def convert_treatmentN_to_treatmentStr(data):
    """
    

    Parameters
    ----------
    data : pd.DataFrame with target_drug column.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    treatment_number = list(np.unique(data.target_drug))
    
    treatments = ["Control",
    "8-zaguanine",
    "Raltitrexed",
    "Topotecan",
    "Floxuridine",
    "Nutlin",
    "Dasatinib",
    "Gefitinib",
    "Vincristine",
    "Bortezomib"] 
    return dict(zip(treatment_number, treatments))
    

def gen_cell_lines_states_replicates():
    """
    """
    cell_lines = ["A549", "MCF7", "RKO"]
    states = ["D", "S"]
    replicates = ["Rep1", "Rep2", "Rep3"]
    var_list = []
    for i in range(len(cell_lines)):
        for j in range(len(states)):
            for k in range(len(replicates)):
                var_str = ""
                var_str += cell_lines[i] + "_"
                var_str += states[j] + "_"
                var_str += replicates[k]
                var_list.append(var_str)
    return var_list



def create_full_data(data_file = "proteinGroups_tryptic_head2000.csv", treshold = 1, proteinID_index = True):
    """
    Creates a matrix with onlye Reporter intensity corrected_treatment_cellLine_replicates
    
    MODIFY HERE TO CHANGE 
        Reporter Intensity
        Peptide tresholding
        
    """
    treshold = treshold
    df = pd.read_csv(data_file, sep = "\t")
    reporter_intensity_corrected = add_prefix_with_treatments()
    unique_peptides = add_prefix()
    
    data = pd.DataFrame()
    
    for i in unique_peptides:
        for j in reporter_intensity_corrected:
            if i.split()[-1] == j.split()[-1]:
                #print(i + "\t - \t" + j)
                data[df[j].name] = treshold_by_peptides(df[j], df[i], treshold = treshold)
    # rename index
    if proteinID_index == True:
        data = data.rename(index=dict(zip(list(range(len(data))), list(df["Protein IDs"]))))
    return data


def log2FC_data(data):
    """
    Converts full data to log2FC values. Use non-normalized values from create_full_data
    """
    log2FC_df = pd.DataFrame()
    for i in range(0,len(data.columns),10):
        i = i
        data_subset = data[data.columns[i:i+10]]
        log_data = data_subset.apply(np.log2)
    
        new_df = pd.DataFrame()
        for j in range(len(log_data.columns)):
            tmp_col = log_data.iloc[:, j].name
            tmp_df = log_data.iloc[:,0] - log_data.iloc[:,j]
            new_df[tmp_col] = tmp_df
            
        log2FC_df = log2FC_df.append(new_df.T)
    log2FC_df = log2FC_df.T
    return log2FC_df

def treshold_by_peptides(intensities, peptide_counts, treshold = 1):
    """
    Function treshold by amount of peptides in samples
    """
    #intensities = df["Reporter intensity corrected 0 A549_D_Rep1"]
    #booleans = (df["Unique peptides A549_D_Rep1"] > treshold)
    
    intensities = intensities
    peptide_counts = peptide_counts
    treshold = treshold
    
    booleans = (peptide_counts > treshold)
    res = intensities * booleans.apply(lambda x: np.nan if x == False else True)
    return res


def add_prefix_with_treatments(prefix = "Reporter intensity corrected", treatments = 10):
    """
    Generate a list with unique peptides and treatments for each unique peptide.
    """
    var_list = gen_cell_lines_states_replicates()
    treatments = [i for i in range(treatments)]
    prefix = prefix
    res_list = []
    for i in var_list:
        for j in treatments:
            unit_str = prefix + " "
            unit_str += str(j) + " "
            unit_str += i
            res_list.append(unit_str)
    return res_list


def add_prefix(prefix = "Peptides"):
    """
    Generates a list of unique peptides and adds a prefix to them.
    """
    var_list = gen_cell_lines_states_replicates()
    prefix = prefix
    res_list = []
    for i in var_list:
        unit_str = prefix + " "
        unit_str += i
        res_list.append(unit_str)
    return res_list

   

def subsetting_data_columns(data, var_name = "Reporter intensity corrected", selected_cell_lines = ["A549", "RKO"], selected_replicates = ["1", "2"], selected_states = ["S"], selected_treatments = ["1", "2", "3"] ):
    """
    Function to generate subset columns of data.
    
    data should be output from create_full_data()
    var_name = <string var name from df>
    selected_cell_lines - <list of cell lines>
    selected_replicates - <list of replicates (str(int))>
    selected_states - <list of states>
    selected_treaments <list of treaments (str(int))>
    """
    data = data
    # split the data headers
    var_name = var_name
    selected_cell_lines = selected_cell_lines
    selected_replicates = selected_replicates
    selected_states = selected_states
    selected_treatments = selected_treatments
    
    selected_subset = []
    for i in data.columns:
        
        cell_line = i.split()[-1].split("_")[0] 
        state = i.split()[-1].split("_")[1]
        replicate = i.split()[-1].split("_")[2][-1]
        treatment = i.split()[-2]
    
        subset_str = var_name + " " 
        if treatment in selected_treatments:
            subset_str += treatment + " "
        else:
            continue
        
        if cell_line in selected_cell_lines:
            subset_str += cell_line + "_"
        else:
            continue
        
        if state in selected_states:
            subset_str += state + "_"
        else:
            continue
        
        if replicate in selected_replicates:
            subset_str += "Rep"+replicate

        else:
            continue
        selected_subset.append(subset_str)
    return selected_subset



def split_data(data):
    """
    Splits the data columns from create_full_data to np.arrays
    """
    data_values = []
    target_drugs = []
    cell_lines = []
    states = []
    replicates = []
    
    for i in data.columns:
        target_drug = i.split()[-2]
        cell_line = i.split()[-1].split("_")[0]
        state = i.split()[-1].split("_")[1]
        replicate = i.split()[-1].split("_")[2]
        
        data_values.append(list(data[i]))
        target_drugs.append(int(target_drug))
        cell_lines.append(cell_line)
        states.append(state)
        replicates.append(replicate)
        
    data_values = np.array(data_values)
    target_drugs = np.array(target_drugs)
    cell_lines = np.array(cell_lines)
    states = np.array(states)
    replicates = np.array(replicates)
    
    return data_values, target_drugs, cell_lines, states, replicates

def create_data_melted_peptide_tryptic_file():
    print("fill in")
    return

protein_file = "proteinGroups tryptic.txt"
peptide_file = "peptides tryptic.txt"


import time 

df_protein = pd.read_csv(protein_file, sep = "\t")
df = pd.read_csv(peptide_file, sep = "\t")

treshold = treshold
df = pd.read_csv(data_file, sep = "\t")


reporter_intensity_corrected = add_prefix_with_treatments(prefix = "Reporter intensity corrected", treatments = 10)
reporter_intensity_count = add_prefix_with_treatments(prefix = "Reporter intensity count", treatments = 10)

cols = ["Leading_razor_protein", "Proteins", "Reporter_intensity_count",
        "Gene_names", "Charges", "Missed_Cleaveges", "PEP", "Score", 
        "Cell_line", "Treatment", "State", "Replicate", "Reporter_intensity_corrected"]
rows = []
t_proc_start = time.time()
for i in range(len(reporter_intensity_corrected)):
    print(str(i) + "/" + str(len(reporter_intensity_corrected)))
    print(time.time() - t_proc_start)
    for j in range(len(df[reporter_intensity_corrected[i]])):
        intensity = df[reporter_intensity_corrected[i]][j]
        count = df[reporter_intensity_count[i]][j]
        PEP = df["PEP"][j]
        score = df["Score"][j]
        charges = df["Charges"][j] # Do we need this?
        missed_cleaveges = df["Missed cleavages"][j]
        proteins = df["Proteins"][j]
        razor_protein = df["Leading razor protein"][j]
        gene = df["Gene names"][j]
        
        split_ = reporter_intensity_corrected[i].split()
        split__ = split_[-1].split("_")
        
        cell_line = split__[0]
        treatment = split_[-2]
        state =  split__[1]
        replicate = split__[2]
        
        row = [razor_protein, proteins, count, gene, charges, missed_cleaveges, PEP, 
               score, cell_line, treatment, state, replicate, intensity]
        rows.append(row)
t_proc_end = time.time()
print(t_proc_end - t_proc_start)


# peptide tresholding and q-value tresholding 
# protein name

data =



    treshold = treshold
    df = pd.read_csv(data_file, sep = "\t")
    reporter_intensity_corrected = add_prefix_with_treatments()
    unique_peptides = add_prefix()
    
    data = pd.DataFrame()
    
    for i in unique_peptides:
        for j in reporter_intensity_corrected:
            if i.split()[-1] == j.split()[-1]:
                #print(i + "\t - \t" + j)
                data[df[j].name] = treshold_by_peptides(df[j], df[i], treshold = treshold)
    # rename index
    if proteinID_index == True:
        data = data.rename(index=dict(zip(list(range(len(data))), list(df["Protein IDs"]))))
    return data




