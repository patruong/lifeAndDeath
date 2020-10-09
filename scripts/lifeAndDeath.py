#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 12:07:58 2020

@author: ptruong
"""


import pandas as pd
import numpy as np
import seaborn as sns

def get_value_multiIdx_df(df, idx_val, value):
    return df.iloc[df.index.get_level_values(idx_val) == value]


data = pd.read_csv("proteinGroups_tryptic_reporterIntensities_noNorm.csv", index_col = 0)
meta = pd.read_csv("proteinGroups_tryptic_sampleMetadata_noNorm.csv", index_col = 0)

data_col = data.columns
tuples = list(zip(data_col.values, meta.treatment.values, meta.cell_line.values, meta.state.values, meta.replicate.values))
index = pd.MultiIndex.from_tuples(tuples, names=['data', 'treatment', "cell_line", "state", "repl"])

df = pd.DataFrame(data.values.T, index = index, columns = data.index)


df_control = df.iloc[df.index.get_level_values("treatment") == "Control"]
df_logFc = pd.DataFrame()
for treatment in meta.treatment.unique():
    df_treatment = df.iloc[df.index.get_level_values("treatment") == treatment]
    #Note imputing -np.inf values with 0
    array_logFc = np.log2(df_control).replace(-np.inf, 0).values - np.log2(df_treatment).replace(-np.inf, 0).values
    idx = df_treatment.index
    col = df_treatment.columns
    df_tmp = pd.DataFrame(array_logFc, index = idx, columns = col) # log2Fc control vs treatment
    df_logFc = df_logFc.append(df_tmp)
    
# Count up vs down regulation for cell treatment fig.4 in paper. Fold-change differential expression.


#############################
# 1. Fold-change differential expression of all alive vs dead cells.
#############################

fold_change_eval = 1 # cut-off treshold

n_larger_than_treshold = (df_logFc > fold_change_eval).T.sum()
n_smaller_than_treshold = (df_logFc < -fold_change_eval).T.sum()

total_count = n_larger_than_treshold + n_smaller_than_treshold    

total_count_dead = total_count.iloc[total_count.index.get_level_values("state") == "D"]
total_count_alive = total_count.iloc[total_count.index.get_level_values("state") == "S"]

total_count_dead.sum()
total_count_dead.std()

total_count_alive.sum()
total_count_alive.std()

print("""
      Dead cells have 141 044 proteins differentially expressed between control and treatments.
      Alive cells have 30 945 proteins differentially expressed between control and treatments.
      """)

df_life_death = pd.DataFrame(np.array([total_count_alive.values, total_count_dead.values]).T, columns=["alive","dead"])
df_life_death.melt(var_name='state', value_name='Measurement')
df_life_death = df_life_death.reset_index().melt(id_vars='index').drop("index", axis = 1).rename({"variable":"state", "value":"count"}, axis = 1)
ax = sns.barplot(x="state", y="count", data=df_life_death)

"""
There seems to be alot more differentiating proteins when we go from 
control group to dead than control group treated and alive. 

There are some cell processes happening when we go from alive to dead, 
so this result is expected. 

We can remove the overlapping proteins to look at what proteins are only
related to the state.
"""

#############################
# 2. Fold-change differential expression of alive vs dead cells for different treatment.
#############################

cell_line_count = get_value_multiIdx_df(total_count, "cell_line", "A549")
cell_line_treatment_count = get_value_multiIdx_df(cell_line_count, "treatment", "Raltitrexed")
cell_line_treatment_state_count = get_value_multiIdx_df(cell_line_treatment_count, "state", "S")
cell_line_treatment_state_count
cell_line_treatment_state_count.sum()

variables = total_count.index.names
states = meta.state.unique()
treatments = meta.treatment.unique()
cell_lines = meta.cell_line.unique()

treatments_ = [] # abbreviations
for i in treatments:
    treatments_.append(i[0])
cell_lines_ = [] # abbreviations
for i in cell_lines:
    cell_lines_.append(i[0])


x_vars = [] # cell_line + treatment abbreviated
x_vars_ = [] # cell_line + treatment non-abbreviated for mapping
x_seps = [] # live and dead
y_vars = [] # counts
for state in states:
    for i in range(len(treatments)):
        treatment = treatments[i]
        for j in range(len(cell_lines)):
            cell_line = cell_lines[j]
            
            cell_line_count = get_value_multiIdx_df(total_count, "cell_line", cell_line)
            cell_line_treatment_count = get_value_multiIdx_df(cell_line_count, "treatment", treatment)
            cell_line_treatment_state_count = get_value_multiIdx_df(cell_line_treatment_count, "state", state).sum()
            
            x_var = cell_lines_[j] + "_" + treatments_[i]
            x_var_ = cell_lines[j] + "_" + treatments[i] # For mapping
            x_sep = state
            y_var = cell_line_treatment_state_count
            
            x_vars.append(x_var)
            x_vars_.append(x_var_) # For mapping
            x_seps.append(x_sep)
            y_vars.append(y_var)


df_treatment_cell_line = pd.DataFrame(np.array([x_vars, x_seps, y_vars]).T, columns=["cell_line_treatment", "state", "n"])
df_treatment_cell_line.n = pd.to_numeric(df_treatment_cell_line.n)

g = sns.catplot(
    data=df_treatment_cell_line, kind="bar",
    x="cell_line_treatment", y="n", hue="state",
    ci="sd", palette="dark", alpha=.6, height=6
)
g.despine(left=True)
g.set_axis_labels("{cell_line}_{treatment}", "count")
g.fig.suptitle('log2Fc differential expression | cutoff ' + "-" + str(fold_change_eval) + " < log2FC(x) < " + str(fold_change_eval))

#############################
# 3. Checking which cell lines and treatment contain  -treshold < log2FC < treshold.
#############################

death_life_ratios = np.array(list((df_treatment_cell_line[df_treatment_cell_line.state == "D"].n))) / np.array(list((df_treatment_cell_line[df_treatment_cell_line.state == "S"].n)))

df_death_life_ratios = pd.DataFrame(np.array([df_treatment_cell_line[df_treatment_cell_line.state == "D"].cell_line_treatment, death_life_ratios]).T, columns = ["cell_line_treatment","n"])
df_death_life_ratios.n = pd.to_numeric(df_death_life_ratios.n)

# toDo Convert abbrevs to full strings so we can read and write about it!
col_mapper = dict(zip(x_vars, x_vars_))
df_death_life_ratios["cell_line_treatment_map"] = df_death_life_ratios.cell_line_treatment.map(col_mapper)


# Ratios of treshold
treshold = 3
df_death_life_ratios[df_death_life_ratios.n > treshold] #dead has larger than treshold ratio
df_death_life_ratios[df_death_life_ratios.n < 1/treshold] # survive has larger than treshold ratio

#############################
# 4. Check if there is a difference between replicates, split surviving and dead cells between replicated
#############################
rep_count = get_value_multiIdx_df(total_count, "repl", "Rep1")
rep_state_count = get_value_multiIdx_df(rep_count, "state", "S")


variables = total_count.index.names
states = meta.state.unique()
replicates = meta.replicate.unique()

x_vars = [] # replicates
x_seps = [] # surviving and death
y_vars = [] # count

for state in states:
    for rep in replicates:
        rep_count = get_value_multiIdx_df(total_count, "repl", rep)
        rep_state_count = get_value_multiIdx_df(total_count, "state", state).sum()
        
        x_var = rep
        x_sep = state
        y_var = rep_state_count
        
        x_vars.append(x_var)
        x_seps.append(x_sep)
        y_vars.append(y_var)

df_rep = pd.DataFrame(np.array([x_vars, x_seps, y_vars]).T, columns = ["replicate", "state", "n"])
df_rep.n = pd.to_numeric(df_rep.n)

g = sns.catplot(
    data=df_rep, kind="bar",
    x="replicate", y="n", hue="state",
    ci="sd", palette="dark", alpha=.6, height=6
)
g.despine(left=True)
g.set_axis_labels("replicate", "count")
g.fig.suptitle('log2Fc differential expression | cutoff ' + "-" + str(fold_change_eval) + " < log2FC(x) < " + str(fold_change_eval))

"""
Replicates show equal amount of larger than logFC treshold and lower than logFC treshold between replicates...
Which is good, because we expect that the replicates behave the same.
"""

#######################
# Fixing the renaming #
#######################


df_logFc.columns[0] 
















    