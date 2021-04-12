#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:48:06 2021

@author: ptruong
"""


import os 
import pandas as pd
import numpy as np

from skbio.stats.composition import closure
from skbio.stats.composition import multiplicative_replacement
from skbio.stats.composition import clr

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt


from get_columns import get_cell_line_state_replicate, get_base_cols_peptide, get_all_reporter_intensity_correct, get_reporter_intensity_without_control
from column_mapper import col_to_treatment_mapper, treatment_nomenclature_map_dict, col_to_cell_line_mapper, col_to_state_mapper, col_to_rep_mapper
pd.options.mode.chained_assignment = None
os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")
df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")


#cell_lines, states, replicates = get_cell_line_state_replicate()
base_cols = get_base_cols_peptide()
reporter_intensity_corrected_cols = get_all_reporter_intensity_correct()

cell_lines, states, replicates = [], [], []
for x in reporter_intensity_corrected_cols:
    cell_line = x.split(" ")[-1].split("_")[0]
    state = x.split(" ")[-1].split("_")[1]
    replicate = x.split(" ")[-1].split("_")[2]
    cell_lines.append(cell_line), states.append(state), replicates.append(replicate)
cell_lines = list(dict.fromkeys(cell_lines).keys())
states = list(dict.fromkeys(states).keys())
replicates = list(dict.fromkeys(replicates).keys())

df_base = df_raw[get_base_cols_peptide()]


def col_to_mIdx(df_int):
    treatment_mapper = lambda x : x.split(" ")[3]
    specie_mapper = lambda x : x.split(" ")[4].split("_")[0]
    state_mapper = lambda x : x.split(" ")[4].split("_")[1]
    replicate_mapper = lambda x : x.split(" ")[4].split("_")[2]
    batch_mapper = lambda x : x.split(" ")[-1]
    experiment_mapper = lambda x: x.split(" ")[3] + "_" + x.split(" ")[-1] 
    sample_mapper = lambda x: x.split(" ")[4].split("_")[0] + "_" + x.split(" ")[4].split("_")[1] + "_" +  x.split(" ")[3]
    
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



def intensities_to_midx_df(df_int):
    """
    df_int = df[reporter_intensity_corrected_cols]

    example on selecting from levels on multiIndex:
        df_int.iloc[:, df_int.columns.get_level_values(level=1) == "S"]
        
    Example used for this build.
    
    Column to midx
    https://riptutorial.com/pandas/example/18696/how-to-change-standard-columns-to-multiindex
    
    Select columns in multiIdx
    https://stackoverflow.com/questions/25189575/pandas-dataframe-select-columns-in-multiindex
    
    # build hierarchical df
    
    df_ = pd.DataFrame(np.random.randn(2,3), columns=['a','b','c'])
    midx = pd.MultiIndex(levels=[['zero', 'one'], ['x','y']], codes=[[1,1,0,],[1,0,1,]])
    df_.columns = midx
    
    
    list1 = ['a','b','c','d'] # values
    list2 = ['d','a', "c", "d"] # value list
    print([list2.index(x) for x in list1 if x in list2])

    """
    midx = col_to_mIdx(df_int)
    df_int.columns = midx
    return df_int

def drop_zero_row(df):
    return df[(df.T != 0).any()]
    
# filter df
df_base.columns
df = df_raw[df_raw.PEP < 0.05] # 5% PEP removed 9578 peptides
df = df[df["Missed cleavages"] < 3] # 0 removed

df_int = df[reporter_intensity_corrected_cols]
df_int = df_int.drop_duplicates() # Removing duplicate rows, these are most likely 0 rows. 2126 rows dropped
df_int = drop_zero_row(df_int) # dropped last zero-row
midx = col_to_mIdx(df_int)

# Raw
df_int = intensities_to_midx_df(df_int)
df_signal = (df_int != 0)
# log2
#df_int = df_int.replace(0, np.nan)
#df_int = np.log2(df_int)


#######################
# BATCH NORMALIZATION #
#######################


df_int.sum() #check intensity sum of each channel

target = df_int.sum().mean()
norm_fac = target/df_int.sum()
df_norm = df_int*norm_fac

df_norm.sum() #check again intensity sum of each channel.

#df_norm_log2 = df_norm.replace(0, np.nan)
#df_norm_log2 = np.log2(df_norm_log2)


batch = pd.DataFrame(df_norm.columns.get_level_values(4)).T.values

from combat.pycombat import pycombat
data_corrected = pycombat(df_norm.fillna(0),batch[0])

#################################
# Generate diffacto input files #
#################################

data_corrected["sequence"] = df_raw.Sequence
data_corrected = data_corrected.set_index("sequence")

df_signal["sequence"] = df_raw.Sequence
df_signal = df_signal.set_index("sequence")

data_corrected_true_signal = data_corrected*df_signal

diffacto_input = pd.DataFrame(data_corrected.values, index = data_corrected.index.values, columns = data_corrected.columns.get_level_values("experiment").values)
diffacto_input.to_csv("peptide_tryptic_combat_diffacto_input.csv", sep = ",")

sample_rep_name = data_corrected.columns.get_level_values("experiment").values
sample_name = data_corrected.columns.get_level_values("sample").values

pd.DataFrame([sample_rep_name, sample_name]).T.to_csv("peptide_tryptic_combat_diffacto_input_sample.lst", index = False, header = False, sep = "\t")
diffacto_input = pd.read_csv("peptide_tryptic_combat_diffacto_input.csv", sep = ",", index_col = 0)

def check_if_row_elements_equal(df):
    return df.eq(df.iloc[:, 0], axis=0)

#####################
# read in diffacto ##
#####################
    
diffacto = pd.read_csv("diffacto_combat_output.protein.txt", sep = "\t")
diffacto_FDR = pd.read_csv("diffacto_combat_output.protein.FDR", sep = "\t")
diffacto["MCFDR"] = diffacto_FDR["MCFDR"]
diffacto = diffacto[diffacto["MCFDR"] < 0.05] #treshold on monte-carlo FDR
diffacto = diffacto[diffacto["S/N"] > 0] # treshold on signal-to-noise ratio


###########
# PCA #####
###########

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import seaborn as sns



def pca_plot_ax(axes, row, col, df_int, classification = "state", marker = "replicate", seperator = "cell_line", seperate_by = "A549",
            title = None):
    # splice df_int for pca now
    #df_pca = df_int[df_int[seperator] == seperate_by]
    df_pca = df_int.iloc[:,df_int.columns.get_level_values(seperator) == seperate_by]
    #df_pca = df_pca.fillna(0)
    df_pca = df_pca.fillna(0).T # transpose because other did not work
    #weights = df_pca.notna()*1
    
    features = df_int.index.values
    
    # Separating out the features
    x = df_pca.loc[:, features].values
    
    # Separating out the target
    # classification = "state # Choose target - treatment, cell_line, state here
    # marker = "rep" #marker
    # y = df_pca.loc[:,[classification]].values
    
    # Standardizing the features
    x = StandardScaler().fit_transform(x)
    
    # Missing value impuration
    imputer = SimpleImputer(missing_values=np.nan, strategy="constant", fill_value = 0)
    x = imputer.fit_transform(x)
    
    # pca
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ["principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                              'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1])])
    #finalDf = pd.concat([principalDf, df_pca[[classification]].reset_index()[classification]], axis = 1)
    finalDf = principalDf
    finalDf[classification] = df_pca.index.get_level_values(classification)
    if marker != False:
        #finalDf = pd.concat([finalDf, df_pca[[marker]].reset_index()[marker]], axis = 1)
        finalDf[marker] = df_pca.index.get_level_values(marker)
    #Visualization
    #fig, ax = plt.subplots(figsize=(16,10))
    if marker != False:
        sns.scatterplot(ax=axes[row,col] , data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]),
                        y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]), 
                        hue = classification, style = marker, s = 100)
    else:
        sns.scatterplot(ax=axes[row,col], data = finalDf, x = "principal component 1 [{:.2f} explained]".format(pca.explained_variance_ratio_[0]), 
                        y = 'principal component 2 [{:.2f} explained]'.format(pca.explained_variance_ratio_[1]), 
                hue = classification, s = 100)
    axes[row, col].set_title(seperate_by)


fig, axes = plt.subplots(3, 3, figsize = (30,24))

cell_lines = ["A549", "RKO", "MCF7"]
color_markers = [["state", "replicate"], ["treatment", "state"], ["treatment", "replicate"]]
for row in range(3):
    color_marker = color_markers[row]
    for col in range(3):
        cell_line = cell_lines[col]
        pca_plot_ax(axes, row, col, df_norm, classification = color_marker[0], marker = color_marker[1], seperator = "cell_line",
                    seperate_by = cell_line, title = None)



### Make the analysis on protein level... with the same type of script

# summarize peptide to protein
            


import time
cell_line = "A549"
state = "D"
treatment = "0"
replicate = "Rep1"

protein = df_int.index.unique()[0]


start = time.time()


def top3_peptide_to_protein_summarization(df_int, cell_line, state, treatment, replicate):
    col_val = "_".join([cell_line, state, replicate])
    col = ("Reporter intensity corrected " + treatment + " " + col_val)
    
    value_array = []
    prot_array = []
    proteins = df_int.index.unique()

    for protein in proteins:
        df_prot = df_int[df_int.index == protein]
        df_val = df_prot[cell_line][state][treatment][replicate].sort_values(by = col_val, ascending = False).head(3)
        if len(df_val) > 2:
            protQuant = df_val.mean().values[0]
        else:
            protQuant = np.nan
    
        prot_array.append(protein)
        value_array.append(protQuant)
    
    df_protSum = pd.DataFrame(value_array, index = prot_array, columns = [col])
    return df_protSum

start = time.time()
iteration = 0
df_all = pd.DataFrame()
df_int["protein"] = df_raw["Leading razor protein"] 
df_int = df_int.set_index("protein")

for cell_line in cell_lines:
    for state in states:
        for treatment in [str(treatment_num) for treatment_num in range(10)]:
            for replicate in replicates:
                print("_".join([str(iteration),cell_line,state,treatment,replicate]))
                iteration += 1
                df_proteinQuantity = top3_peptide_to_protein_summarization(df_int, cell_line, state, treatment, replicate)
                df_all = pd.concat([df_all, df_proteinQuantity], axis = 1)
                print(time.time()-start)
df_all.to_csv("df_all.csv", sep = "\t")                
end = time.time()
print(end-start)                
               



################################
# protein sums normalization ###
################################


df_all = pd.read_csv("protein_summarized.csv", sep = "\t", index_col = 0)

# Raw
df_all = intensities_to_midx_df(df_all)


df_all.sum() #check intensity sum of each channel

target = df_all.sum().mean()
norm_fac = target/df_all.sum()
df_norm_prot = df_all*norm_fac

df_norm_prot.sum() #check again intensity sum of each channel.

batch = pd.DataFrame(df_norm_prot.columns.get_level_values(4)).T.values

from combat.pycombat import pycombat
#data_corrected = pycombat(df_norm_prot,batch)
data_corrected = pycombat(df_norm.fillna(0),batch[0])







##################### THiNK ABOUT THIS... maybe aitchison before ComBat?
########################################
# Aitchison multiplicative_replacement #
########################################
data_corrected.sum()
df_int.sum()
df_aitchison = multiplicative_replacement(df_int)
df_aitchison = pd.DataFrame(df_int, columns = midx)


def aitchison_transform(df):
    """
    Aitchison tranformation on df.
    
    df should consist of all samples tagged together in one channel (i.e. A549_S_rep1 etc.)
    """
    df_aitchison = multiplicative_replacement(df)
    #df_aitchison = closure(df)
    df_idx = df.index
    df_col = df.columns
    df_aitchison = pd.DataFrame(df_aitchison, index = df_idx, columns = df_col)
    return df_aitchison

df_aitchison = pd.DataFrame()

start = time.time()
for cell_line in cell_lines:
    for state in states:
        for rep in replicates:
            df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == ("_".join([cell_line, state, rep]))]
            df_part = df_part[(df_part.T != 0).any()]
            df_part = aitchison_transform(df_part)
            df_aitchison = pd.concat([df_aitchison, df_part], axis = 1 )
            print(time.time()-start)
end=time.time()
print(end-start)


cell_line = "A549"
state = "S"
rep = "Rep1"
df1 = df_int.iloc[:, df_int.columns.get_level_values(4) == ("_".join([cell_line, state, rep]))]
df1 = df1[(df1.T != 0).any()]
df_aitchison = multiplicative_replacement(df1)
df_idx = df1.index
df_col = df1.columns
df_aitchison = pd.DataFrame(df_aitchison, index = df_idx, columns = df_col)

df1 = df_int.iloc[:, df_int.columns.get_level_values(4) == ("_".join([cell_line, state, rep]))]
df1 = df1[(df1.T != 0).any()]

df2 = df_int.iloc[:, df_int.columns.get_level_values(4) == ("_".join([cell_line, "D", rep]))]
df2 = df2[(df2.T != 0).any()]

pd.concat([df1, df2], axis = 1 )


a = pd.DataFrame(np.array([[1,2,3], [4,0,2], [4,5,6], [3,3,3]]))
b = closure(a)
b = multiplicative_replacement(a)
pd.DataFrame(b)
pd.DataFrame(b).sum(axis=1)

###########################################################
# PLOT histogram - investigate error in A549_D aitchison ##
###########################################################
import seaborn as sns

df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                      ("_".join(["A549", "D", "Rep2"]))]
for i in range(np.shape(df_part)[1]):
    #plt.hist(np.log2((df_part.iloc[:,i]).replace(0,np.nan)), bins = 1000, 
    #         histtype="step", label = df_part.iloc[:,i].name[-1])
    sns.distplot(np.log2((df_part.iloc[:,i]).replace(0,np.nan)), bins = 1000,
                 label = df_part.iloc[:,i].name[-1], kde = True, hist = False)
plt.legend()

ax = np.log2(df_part.replace(0,np.nan)).plot.hist(bins=1000, alpha = 0.5)

df_part = df_part[(df_part.T != 0).any()]
df_ait_part = aitchison_transform(df_part)


for i in range(np.shape(df_ait_part)[1]):
#    plt.hist((df_ait_part.iloc[:,i]).replace(0,np.nan), bins = 1000, 
#             histtype="step", label = df_part.iloc[:,i].name[-1])
    sns.distplot(df_ait_part.iloc[:,i], bins = 1000,
             label = df_part.iloc[:,i].name[-1], kde = True, hist = False)
plt.legend()



###########3###################################
# Plot all aitchison withing batch histogram ##
###############################################

df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                      ("_".join(["A549", "D", "Rep2"]))]
df_part = df_part[(df_part.T != 0).any()]
df_ait_part = aitchison_transform(df_part)

sns.distplot(np.log2(df_part.replace(0,np.nan)))

sns.distplot(np.log2(df_part.replace(0,np.nan)).values.flatten())

# raw - all intensities within batch together
for cell_line in cell_lines:
    for state in states:
        for replicate in replicates:
            df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                      ("_".join([cell_line, state, replicate]))]
            sns.distplot(np.log2(df_part.replace(0,np.nan)),
                         label = cell_line + "_" + state + "_" + replicate,
                         hist = False, kde = True)
plt.legend()
            

# aitchison - all intensities within batch togethre
for cell_line in cell_lines:
    for state in states:
        for replicate in replicates:
            df_part = df_int.iloc[:, df_int.columns.get_level_values("batch") == 
                      ("_".join([cell_line, state, replicate]))]
            df_part = df_part[(df_part.T != 0).any()]
            df_ait_part = aitchison_transform(df_part)
            sns.distplot(df_ait_part,
                         label = cell_line + "_" + state + "_" + replicate,
                         hist = False, kde = True)
plt.legend()
            



sns.distplot(df_ait_part)
sns.distplot(df_ait_part.values.flatten())
########
# PLOT #
########

def plot_intensity_histogram(df_int, min_x = 0, max_x = 30, step = 0.1, title = "title"):
    """
    df_int with midx
    
    plots the step histogram for {cell_line}_{state}
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    col_i = 0
    for cell_line in cell_lines:
        for state in states:
            df_state = df_int.iloc[:, df_int.columns.get_level_values(level=1) == state]
            df_cell_state = df_state.iloc[:, df_state.columns.get_level_values(level=0) == cell_line]
            plot_label = True
            for i in range(np.shape(df_cell_state)[1]):
                if plot_label == True:
                    plt.hist(df_cell_state.iloc[:,i] , bins=np.arange(min_x,max_x, step), 
                             histtype="step", color = colors[col_i], alpha = 0.4,
                             label = cell_line + "_" + state)
                    plot_label = False
                else:
                    plt.hist(df_cell_state.iloc[:,i] , bins=np.arange(min_x,max_x, step), 
                             histtype="step", color = colors[col_i], alpha = 0.4)
            col_i+=1
    plt.title(title)
    plt.legend()
    

def plot_intensity_boxplot(df_int, title = "title"):
    """
    Boxplot details.
    # https://stackoverflow.com/questions/41997493/python-matplotlib-boxplot-color

    df_int with midx
    
    plots the step histogram for {cell_line}_{state}
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    col_i = 0
    pos_i = 0
    for cell_line in cell_lines:
        for state in states:
            df_state = df_int.iloc[:, df_int.columns.get_level_values(level=1) == state]
            df_cell_state = df_state.iloc[:, df_state.columns.get_level_values(level=0) == cell_line]
            plot_label = True
            for i in range(np.shape(df_cell_state)[1]):
                if plot_label == True:
                    plt.boxplot(df_cell_state.iloc[:,i].dropna(), positions = [pos_i], patch_artist=True, 
                                notch=True, medianprops=dict(color="grey"),
                                boxprops=dict(facecolor = colors[col_i]),)
                                # No label yet, so if-else does not matter.......
                    plot_label = False
                else:
                    plt.boxplot(df_cell_state.iloc[:,i].dropna(), positions = [pos_i], patch_artist=True, 
                                notch=True, medianprops=dict(color="grey"),
                                boxprops=dict(facecolor = colors[col_i]))
                pos_i+=1
            col_i+=1
    plt.title(title)
    plt.legend()  

plot_intensity_histogram(np.log2(df_int.replace(0,np.nan)), min_x = 0, max_x = 25, step = 0.1, title = "raw")
plot_intensity_boxplot(np.log2(df_int.replace(0,np.nan)), title = "raw")

plot_intensity_histogram(df_aitchison, min_x = 0, max_x = 1, step = 1/1000, title = "raw")
plot_intensity_boxplot(df_aitchison, title = "raw")

