#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:16:54 2020

@author: ptruong
"""

import pandas as pd
import numpy as np
from helper import *



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


def count_nan(df):
    count_nan = len(df) - df.count()
    return count_nan


def create_full_data(data_file = "proteinGroups_tryptic_head2000.csv", treshold = 1):
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
                print(i + "\t - \t" + j)
                data[df[j].name] = treshold_by_peptides(df[j], df[i], treshold = treshold)
    return data


   

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
 
def x_y_split_data(data):
    """
    Splits data into X and y np.array
    """
    X = []
    y = []
    for i in data.columns:
        target_drug = i.split()[-2]
        y.append(int(target_drug))
        X.append(list(data[i]))
    X = np.array(X)
    y = np.array(y)
    return X, y


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

def normalize(df):
    """
    mean normalization column-wise
    """
    normalized_df=(df-df.mean())/df.std()
    return normalized_df
X, y = x_y_split_data(data)
    

if __name__ == "__main__":    
    filename = "proteinGroups tryptic.txt"
    data = create_full_data(data_file = filename, treshold = 1)
    #subset_cols = subsetting_data_columns(data, var_name = "Reporter intensity corrected", selected_cell_lines = ["A549", "RKO", "MCF7"], selected_replicates = ["1", "2", "3"], selected_states = ["S"], selected_treatments = ["1", "2"] )    
    #subset = data[subset_cols]
    #subset = data
    #X,y = x_y_split_data(subset)

    #subset = subset.fillna(0) # Imput missing value to zero
    #subset = normalize(subset)
    #subset = subset.T
    #subset['label'] = y

    data = data.fillna(0) # Input missing values to zero
    data = normalize(data)
    data_values, target_drugs, cell_lines, states, replicates = split_data(data)
    data = data.T
    data['target_drug'] = target_drugs
    data['cell_line'] = cell_lines
    data['states'] = states
    data['replicates'] = replicates
    
    styles = []
    for i in cell_lines:
        for j in states:
            style = i + "_" + j 
            styles.append(style)
    styles = np.array(styles)
    data["styles"] = styles
    #subset = subset.values
    

# ToDo split data on cell line, treatment and 

###########
# PCA #####
###########

np.random.seed(42)
rndperm = np.random.permutation(data.shape[0])
data = data.reset_index().drop("index", axis = 1)

pca = PCA(n_components=3)
pca_result = pca.fit_transform(subset)
data['pca-one'] = pca_result[:,0]
data['pca-two'] = pca_result[:,1] 
data['pca-three'] = pca_result[:,2]
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))




markers = {"Lunch": "s", "Dinner": "X"}

ax = sns.scatterplot(x="total_bill", y="tip", style="time",
                     markers=markers,
                     scatter_kws = {'facecolors':'none'},
                     data=tips)


# Create marker dict
marker_dict = dict()
fillstyles = []
markers = ["o","1",
           "^","x",
           "s","+"]
count = 0
for i in np.unique(cell_lines):
    for j in np.unique(states):
        marker_for = i+"_"+j        
        marker_dict.update({marker_for: markers[count]})
        count += 1

data["marker"] = data["cell_line"]+"_"+data["states"]




sns.scatterplot(
    x="pca-one", y="pca-two",
    hue="target_drug",
    palette=sns.color_palette("hls", len(np.unique(labels))),
    data=data.loc[rndperm,:],
    legend="full",
    alpha=0.3
)





# PCA plot 2D
plt.figure(figsize=(16,10))
#markers = {'A549':"s", 'MCF7':"X", 'RKO':"o"}
ax = sns.scatterplot(x="pca-one", y="pca-two", 
                     style="marker",
                     hue = "target_drug",
                     palette=sns.color_palette("hls", len(np.unique(labels))),
                     legend="full",
                     alpha=0.7,
                     #markers=markers,
                     data=data.loc[rndperm,:])




ax = plt.figure(figsize=(16,10)).gca(projection='3d')
ax.scatter(
    xs=data.loc[rndperm,:]["pca-one"], 
    ys=data.loc[rndperm,:]["pca-two"], 
    zs=data.loc[rndperm,:]["pca-three"], 
    c=data.loc[rndperm,:]["target_drug"], 
    cmap='tab10',
)
ax.set_xlabel('pca-one')
ax.set_ylabel('pca-two')
ax.set_zlabel('pca-three')
plt.show()





########################
#######################

pca = PCA(n_components=3)
pca_result = pca.fit_transform(subset)
subset['pca-one'] = pca_result[:,0]
subset['pca-two'] = pca_result[:,1] 
subset['pca-three'] = pca_result[:,2]
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

plt.figure(figsize=(16,10))
sns.scatterplot(
    x="pca-one", y="pca-two",
    hue="label",
    palette=sns.color_palette("hls", len(np.unique(labels))),
    data=subset.loc[rndperm,:],
    legend="full",
    alpha=0.3
)

ax = plt.figure(figsize=(16,10)).gca(projection='3d')
ax.scatter(
    xs=subset.loc[rndperm,:]["pca-one"], 
    ys=subset.loc[rndperm,:]["pca-two"], 
    zs=subset.loc[rndperm,:]["pca-three"], 
    c=subset.loc[rndperm,:]["label"], 
    cmap='tab10'
)
ax.set_xlabel('pca-one')
ax.set_ylabel('pca-two')
ax.set_zlabel('pca-three')
plt.show()


data_subset = df_subset[feat_cols].values




