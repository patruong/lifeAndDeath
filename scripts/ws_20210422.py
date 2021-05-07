#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:23:04 2021

@author: ptruong
"""

import os 
os.chdir("/home/ptruong/git/lifeAndDeath/scripts")
from q_value import qvalues

import pandas as pd
import numpy as np
import time

os.chdir("/home/ptruong/git/lifeAndDeath/data/amirata")

df_raw = pd.read_csv("peptides tryptic.csv", sep = "\t")



cell_lines = [ "A549", "MCF7", "RKO"]
states = ["S", "D"]
treatments = [i for i in range(10)]


# Check so that sequence only maps to one protein.
df_raw.groupby("Sequence")["Leading razor protein"].apply(len).mean()

cell_line = "A549"
treatment = 2
state = "S"

#from scipy.stats import ttest_ind

start = time.time()

sample = []
dfs = []
for cell_line in cell_lines:
    for state in states:
        for treatment in treatments:
            #df_log2FC_protQuant = pd.DataFrame()
            protein_array = []
            log2FC_array = []
            
            for protein, peptides in df_raw.groupby("Leading razor protein"):
                score = peptides["Score"]
                reduced = peptides.copy().filter(regex=(f"^Reporter intensity corrected [0{treatment}] {cell_line}_{state}.*"))
                reduced["score"] = score
                #reduced['sum'] = reduced.sum(axis=1)
                ctrl_col = [col for col in reduced if col.startswith("Reporter intensity corrected 0")]
                case_col = [col for col in reduced if col.startswith(f"Reporter intensity corrected {treatment}")]
                
                #reduced.sort_values(by="sum", ascending=False, inplace=True)
                reduced.sort_values(by="score", ascending=False, inplace=True)
                proteinQuant = reduced.head(1).drop("score", axis = 1)
                proteinQuant.rename({proteinQuant.index[0]:protein}, inplace=True)
                
                log_proteinQuant_ctrl = np.log2(proteinQuant.iloc[0:3,:][ctrl_col].replace(0,np.nan))
                log_proteinQuant_case = np.log2(proteinQuant.iloc[0:3,:][case_col].replace(0,np.nan))
                log2FC = log_proteinQuant_ctrl.values - log_proteinQuant_case.values
                
                protein_array.append(protein)
                log2FC_array.append(log2FC[0])
            
            end = time.time()
            print(end-start)
    
            df_protQuant = pd.DataFrame(log2FC_array, columns = case_col, index = protein_array)
            dfs.append(df_protQuant)
            sample.append(cell_line + "_" + state + "_" + str(treatment))

res = pd.concat(dfs, axis = 1)
res.to_csv("ws_20210422.csv", sep = "\t")
res = pd.read_csv("ws_20210422.csv", sep = "\t", index_col = 0)



# Check assumtions of ANOVA...
"""
i) Jag tycker också att du bör börja med att känna på datat. 
    Tex skulle du för en en given cell linje och behandling kunna göra en plott med 
    
x-axel: FC mellan behandlad/obehandlad i levande
y-axel: FC mellan behandlad/obehandlad i döda

ii) Jag tror att du bör ta reda på om ditt data, efter det att du normalisrat det 
    genom att dela med control, verkligen är normalfördelat. Jag misstänker att det 
    inte är det, och att det kanske är så att du inte kan använda en anova alls.

iii) Även om residualerna skulle kunna aproximeras som normalfördelade, så tror 
    jag att du har på tok för mycket batcheffekt för att kunna lita på p-värdena 
    från en Annova annat än innom varje cell-linje och död/levande tillstånd. 
    Dvs kanse att 2 funkar, men förmodligen inte de andra.
    
    
i) För att ta redan på om det är homoskedastiskt
ii) Om det inte är normalfördelat, är det inte därför man använder Kruskal-Wallis 
    Anova? Jag är inte helt säker på vilka konsekvenser Kruskal-Wallin ANOVA har dock.. 
    ska läsa på lite om detta.
iii) Jag trodde att det här med batcheffekterna inte skulle vara något problem 
    efter FC eftersom vi bara kollar på hur många FC det skiljer sig mellan 
    kontrol och behandling inom en kanal. Eller är det så att batch-effekter 
    som dosering av t.ex. behandling etc. innan som kan påverka även FC mellan 
    kontrol och behandling?  
"""

cell_line = "MCF7"
treatment = "1"
S = res.copy().filter(regex=(f"^Reporter intensity corrected [{treatment}] {cell_line}_S.*"))
D = res.copy().filter(regex=(f"^Reporter intensity corrected [{treatment}] {cell_line}_D.*"))
col_mapper = {}

x = "Reporter intensity corrected 1 A549_D_Rep1"

col_mapper = lambda x: x[-4:]
S = S.rename(columns=col_mapper)
D = D.rename(columns=col_mapper)
SD = pd.concat([S.stack() ,D.stack()], axis = 1).rename(columns={0:"S", 1:"D"})

import matplotlib.pyplot as plt

import seaborn as sns
sns.scatterplot(data = SD, x = "S", y = "D")
plt.title("FC S vs D (MCF7 C vs T1)")

# same thing but for summed protein, perform t-test to filter away trash
def shapiro(df):
    shapiro = pd.DataFrame(df.apply(stats.shapiro, axis = 1).tolist(), index = df.index)
    shapiro["q"] = qvalues(shapiro, pcolname="pvalue") 
    return shapiro

shapiro_S = shapiro(S) #Shapiro p > 0.05 indicates normality
shapiro_D = shapiro(D)
shapiro_S.q.min() 
shapiro_D.q.min()

ttest = pd.DataFrame(stats.ttest_ind(S,D, axis = 1), columns = S.index, index = ["t", "p"]).T
ttest["q"] = qvalues(ttest, pcolname="p")
ttest["p"].hist()
plt.title("t-test, (RKO C vs T3)")

fig, ax = plt.subplots()
S.stack().hist(bins = 100)
plt.title("S all FC (RKO C vs T3)")

fig, ax = plt.subplots()
D.stack().hist(bins = 100)  
plt.title("D all FC (RKO C vs T3)")
D.hist(bins=100)

fig, ax = plt.subplots()
SD.plot(kind = "hist", bins = 100, grid = True, alpha = 0.7)
plt.title("All FC (MCF7 C vs T1)")

SD_mean = pd.concat([S.mean(axis=1), D.mean(axis=1)], axis = 1).rename(columns={0:"S", 1:"D"})

sns.scatterplot(data = SD_mean, x = "S", y = "D")
plt.title("FC S vs D, avg. protQuant")

import scipy.stats as stats


### Make subplot with all samples...ö
### Make matrix with full shapiro test
### Check ANOVA assumptions
from statsmodels.stats.outliers_influence import variance_inflation_factor

X = S.dropna()

variance_inflation_factor(X.values, 2)
# some proteins will always be up-regulated on dead and vice versa.

# look at the same thing for summed proteins.

#####
cell_line = "MCF7"
treatment = "1"
S = res.copy().filter(regex=(f"^.*_S.*"))
D = res.copy().filter(regex=(f"^.*_D.*"))

col_mapper = lambda x : x.split(" ")[-2] + "_" + x[-4:] + "_" + x.split(" ")[-1].split("_")[0]

S = S.rename(columns=col_mapper)
D = D.rename(columns=col_mapper)
SD = pd.concat([S.stack() ,D.stack()], axis = 1).rename(columns={0:"S", 1:"D"})

sns.scatterplot(data = SD, x = "S", y = "D")
plt.title("FC S vs D (MCF7 C vs T1)")



shapiro_S = shapiro(S) #Shapiro p > 0.05 indicates normality
shapiro_D = shapiro(D)
shapiro_S.q.min() 
shapiro_D.q.min()



 







