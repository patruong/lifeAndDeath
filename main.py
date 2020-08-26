#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:37:51 2020

@author: ptruong
"""

from __future__ import print_function
import time
import numpy as np
import pandas as pd
from sklearn.datasets import fetch_mldata
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

from preprocessor import *
from stats import *
from plotting import *

if __name__ == "__main__":    
    filename = "proteinGroups tryptic.txt"
    data = create_full_data(data_file = filename, treshold = 1)

    data = data.fillna(0) # Input missing values to zero
    data = log2FC_data(data)
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.fillna(0) # Input missing values to zero
    #data = normalize(data)
    data_values, target_drugs, cell_lines, states, replicates = split_data(data)
    data_vals = data.T # create a data value df
    data['target_drug'] = target_drugs
    data = data.T
    data['cell_line'] = cell_lines
    data['states'] = states
    data['replicates'] = replicates
    data["marker"] = data["cell_line"]+"_"+data["states"]
    target_drug_map = convert_treatmentN_to_treatmentStr(data)
    data["treatment"] = data.replace({"target_drug": target_drug_map}).target_drug

###########
# PCA #####
###########

np.random.seed(42)
rndperm = np.random.permutation(data.shape[0])
data_vals = data_vals.reset_index().drop("index", axis = 1)
data = data.reset_index().drop("index", axis = 1)


pca = PCA(n_components=3)
pca_result = pca.fit_transform(data_vals)
data['pca-one'] = pca_result[:,0]
data['pca-two'] = pca_result[:,1] 
data['pca-three'] = pca_result[:,2]
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))


# PCA plot 2D
plt.figure(figsize=(16,10))
#markers = {'A549':"s", 'MCF7':"X", 'RKO':"o"}
ax = sns.scatterplot(x="pca-one", y="pca-two", 
                     style="marker",
                     hue = "treatment",
                     palette=sns.color_palette("hls", len(np.unique(data.target_drug))),
                     legend="full",
                     alpha=0.7,
                     #markers=markers,
                     data=data.loc[rndperm,:])

# PCA plot 3D
ax = plt.figure(figsize=(16,10)).gca(projection='3d')
ax.scatter(
    xs=data.loc[rndperm,:]["pca-one"], 
    ys=data.loc[rndperm,:]["pca-two"], 
    zs=data.loc[rndperm,:]["pca-three"], 
    c=data.loc[rndperm,:]["treatment"], 
    cmap='tab10',
)
ax.set_xlabel('pca-one')
ax.set_ylabel('pca-two')
ax.set_zlabel('pca-three')
plt.show()

###########
# t-SNE ###
###########

time_start = time.time()
tsne = TSNE(n_components=2, verbose=1, perplexity=100, n_iter=5000)
#tsne = TSNE(n_components=2, verbose=1, perplexity=100, n_iter=300)
tsne_results = tsne.fit_transform(data_vals)
print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start))

data['tsne-2d-one'] = tsne_results[:,0]
data['tsne-2d-two'] = tsne_results[:,1]

# t-SNE plot 2D
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    style="marker",
    hue="target_drug",
    palette=sns.color_palette("hls", len(np.unique(data.target_drug))),
    data=data,
    legend="full",
    alpha=0.7
)

# t-SNE plot 3D <--------- FAILED
time_start = time.time()
tsne = TSNE(n_components=3, verbose=1, perplexity=1000, n_iter=1000)
tsne_results = tsne.fit_transform(data_vals)
print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start))

data['tsne-3d-one'] = tsne_results[:,0]
data['tsne-3d-two'] = tsne_results[:,1]
data['tsne-3d-three'] = tsne_results[:,2]

ax = plt.figure(figsize=(16,10)).gca(projection='3d')
ax.scatter(
    xs=data.loc[rndperm,:]["tsne-3d-one"], 
    ys=data.loc[rndperm,:]["tsne-3d-two"], 
    zs=data.loc[rndperm,:]["tsne-3d-three"], 
    c=data.loc[rndperm,:]["treatment"], 
    cmap='tab10',
)
ax.set_xlabel('pca-one')
ax.set_ylabel('pca-two')
ax.set_zlabel('pca-three')
plt.show()


########
# t-SNE on PCA (dim reduced)






















