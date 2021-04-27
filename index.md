# About

## Relevant Papers
[ProTargetMiner - Saei et al. 2019](https://www.nature.com/articles/s41467-019-13582-8)
[Comparative Proteomics - Saei et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29572246/)



# Blog

## 2021-04-27 Checked anova assumptions.

We need to remove A549 Rep1. It messes stuff up...

## 2021-04-22 top3 or score for protein quant.

Score seems to be more reasonable, because if we take top 3, we could get pepeptides which are 0 (nan) on some columns channels and filter away these peptides. I guess the best scoring peptide should be the most accurate representation of protein quantity.

Why can we take one peptide as protein quantity? Because we are not looking at absolut protein quantities, but in the end we want to perform DE anyway, then we just need to have the same relational values.


## 2021-04-21 ANOVA and Tukey

Tukey HSD discussion....

https://www.researchgate.net/post/Three-means-comparison-by-t-test-or-ANOVA


## 2021-04-20 Top3 again...
https://pubs.acs.org/doi/10.1021/pr800649s

Top3 is top3 peptides per protein in each injection, excluding shared and modified peptides.

MaxQuant info
http://www.coxdocs.org/doku.php?id=maxquant:table:peptidetable

General data info
https://sites.psu.edu/msproteomics/tag/protein-group/

Looking into MANOVA

https://statistics.laerd.com/spss-tutorials/one-way-manova-using-spss-statistics.php

Here are some assumptions for MANOVA.

WHY ANOVA?
https://statistics.laerd.com/statistical-guides/one-way-anova-statistical-guide.php

## 2021-04-19 Top3 and protein summarization error.

Top3 is not the way I do it.

index should not be leading razor protein.

The average should be between technical replicates.


## 2021-04-18 Looking into Top3 and writing minimal script.

Writing minimal script ws_20210418_top3_minimal_example.py

Change top3 computations from >2 to >=2 peptides. The >= is what is used in LFQ Bench in Tenzer.

## 2021-04-14-15 Looked at a lot of histrograms, PCA and restructured the code and wrote report.

As title mentions... A report in google drive..


## 2021-04-13 Protein summarization seem to yield good results.

Pipeline is 

normalize, the ratios of values between samples should be same (pwilmart), to remove batch effect.

Perform aitchison to allow to work with the data in euclidean space.

TOP3 (average of top 3 peptide intensities as protein quantity) to get protein intensities.

raw -> sl -> irs -> tmm -> aitchison -> Top3 

Protein existing in less than 2 samples are discarded.

Going to make the DE protein experiment as well.

## 2021-04-10-12 More histograms, pca and protein summarization.

Kept trying to fix the batch effect... Combat seems to be doing some wierd fill-ins of 0 (NaN) values, that I don't want. They do not make sense. (I should read combat to understand what to make of this).

Recoded the the raw->SL->IRS->TMM normalization and now it seems to work. I think i got the IRS part wrong ealier. 

It works even better with aitchison.

Protein summarization implemented with average of top 3 intensities, and if there are less than 2 peptide intensities, it is discarded. 

alternative method for protein quantification is given by iPQF (for labelled quantification) [Fischer et al. 2016](https://academic.oup.com/bioinformatics/article/32/7/1040/1743612).

Running triqler for FC 1.05 to 2.00 for the dataset. 

I should to pq-plot for triqler and osw data.



## 2021-04-08 Aitchison and IRS Normalization remake in R.

Realized that I might have done the historgrams wrong yesterday.

Thought about how to perform the transformations.

ATM - Aitchison -> SL Norm -> ComBat -> Diffacto

Aitchison seem to center the values.

SL Norm does not to seem to make sense. I will replicate pwilmarts R script to see what is going on.
## 2021-04-06-07 Protein summarization and Diffacto

Wrote script for summarizing proteins based on top3. It's really slow.

Got Diffacto working on the data. Unclear stuff:
- In which order should we normalize, filter and run diffacto.

Run of diffacto without combat yields wierd histograms.

Run of diffacto on ComBat normalized data gives much better histograms. Perhaps because the guassian mixtures get "correct" without the batch... otherwise we make the guassian mixture model on the batches rather than peptide and protein concentration, because we group the samples together and there is a batch effect between replicates.


## 2021-03-28 IRS code
The exact IRS code with normalizing for each batch does not work. Myu old code using geometric average for all batches works better...

I will try comBat tomorrow.

## 2021-03-21-25 Data formatting and batch normalization.

Performed this https://pwilmart.github.io/IRS_normalization/understanding_IRS.html and plotted.



## 2021-03-18-21 Analysis scripts made for proteinGroup tryptic.

jupyter notebook for analysis of proteinGroup.csv is created.

Analysis performed:
- PCA on data seperated on cell lines, states and cell lines & states.
- Batch effect investigation using intensity histograms of data seperated on cell lines & states.

Tested modules for wPCA.
- statsmodels can weight PCA but requires fillna before with ruins the purpose.
- wpca module for python is insanely slow.

I guess weighted PCA with 0 for nan and 1 for non-nan should be same as fillna(0), so I used fillna(0). fillna(0) compares to non fillna(0) make the PCA plot groups "tighter".
ToDo:
- Perform same analysis on peptide data.
- Perform Aitchison transformation and log2FC transformation the "correct" way on this data (before I performed it on bunched together cell lines).


## 2021-03-17 Script for protein analysis

Made script for batching up the proteinGroup file.

Made functions for tresholding for peptide count and sequence coverage.

## 2021-03-15 Batch-effects.

plotted stacked histogrames of the df subsets to see if there was batch effects. It does not seem to be so...

Lukas seem to see that there are signifikant batch effect.

## 2021-03-14 Understanding and wrangling the data.
[MaxQuant - i.e. our file - format info](http://proteomics.fiocruz.br/software/sim/supplementaryfiles/AndromedaResultsPFU/tables.pdf)

score - Andromeda score - how should we use this ?

Problems arising...

- the seperated state has all zero-rows. (probabily, because there is some peptide that is found in either dead or suriving state that is not prevalent in the other.)

When merging some peptide are lost... why?

len(df) = 225715
len(df_z) = 225715
len(df_merged) = 222935

difference is 2781 peptides. Do these matter?

When adding df_base (non reporter_intensity_corrected cols) to the merged. The rows are again 225715. I guess some of the row values are removed in the normalization and clr process? (Does it matter?)

Data merging with aitchison transformation done.

ToDo:

- Do clustering with and without aitchison transformed data.
- Make code for merging data without aitchison transformation -> make normalization. 

## 2021-03-03 Table of experiment
Sent to confirm experiment with Amirata.

| TMT tag | surviving_replicate1                       | surviving_replicate2                       | surviving_replicate3                       | dead_replicate1                            | dead_replicate2                            | dead_replicate3                            |
|---------|--------------------------------------------|--------------------------------------------|--------------------------------------------|--------------------------------------------|--------------------------------------------|--------------------------------------------|
| 126     | Control_1 (Attached cells)                 |  Control_2 (Attached cells)                | Control_3 (Attached cells)                 | Control_1 (Detached cells)                 | Control_2 (Detached cells)                 | Control_3 (Detached cells)                 |
| 127N    | 8-azaguanine on Control_1 (Attached cells) | 8-azaguanine on Control_2 (Attached cells) | 8-azaguanine on Control_3 (Attached cells) | 8-azaguanine on Control_1 (Detached cells) | 8-azaguanine on Control_2 (Detached cells) | 8-azaguanine on Control_3 (Detached cells) |
| 127C    | Raltitrexed on Control_1 (Attached cells)  | Raltitrexed on Control_2 (Attached cells)  | Raltitrexed on Control_3 (Attached cells)  | Raltitrexed on Control_1 (Detached cells)  | Raltitrexed on Control_2 (Detached cells)  | Raltitrexed on Control_3 (Detached cells)  |
| 128N    | Topotecan on Control_1 (Attached cells)    | Topotecan on Control_2 (Attached cells)    | Topotecan on Control_3 (Attached cells)    | Topotecan on Control_1 (Detached cells)    | Topotecan on Control_2 (Detached cells)    | Topotecan on Control_3 (Detached cells)    |
| 128C    | Floxuridine on Control_1 (Attached cells)  | Floxuridine on Control_2 (Attached cells)  | Floxuridine on Control_3 (Attached cells)  | Floxuridine on Control_1 (Detached cells)  | Floxuridine on Control_2 (Detached cells)  | Floxuridine on Control_3 (Detached cells)  |
| 129N    | Nutlin on Control_1 (Attached cells)       | Nutlin on Control_2 (Attached cells)       | Nutlin on Control_3 (Attached cells)       | Nutlin on Control_1 (Detached cells)       | Nutlin on Control_2 (Detached cells)       | Nutlin on Control_3 (Detached cells)       |
| 129C    | Dasatanib on Control_1 (Attached cells)    | Dasatanib on Control_2 (Attached cells)    | Dasatanib on Control_3 (Attached cells)    | Dasatanib on Control_1 (Detached cells)    | Dasatanib on Control_2 (Detached cells)    | Dasatanib on Control_3 (Detached cells)    |
| 130N    | Gefitinib on Control_1 (Attached cells)    | Gefitinib on Control_2 (Attached cells)    | Gefitinib on Control_3 (Attached cells)    | Gefitinib on Control_1 (Detached cells)    | Gefitinib on Control_2 (Detached cells)    | Gefitinib on Control_3 (Detached cells)    |
| 130C    | Vincristine on Control_1 (Attached cells)  | Vincristine on Control_2 (Attached cells)  | Vincristine on Control_3 (Attached cells)  | Vincristine on Control_1 (Detached cells)  | Vincristine on Control_2 (Detached cells)  | Vincristine on Control_3 (Detached cells)  |
| 131     | Bortezomib on Control_1 (Attached cells)   | Bortezomib on Control_2 (Attached cells)   | Bortezomib on Control_3 (Attached cells)   | Bortezomib on Control_1 (Detached cells)   | Bortezomib on Control_2 (Detached cells)   | Bortezomib on Control_3 (Detached cells)   |

| TMT tag | surviving_replicate1 | surviving_replicate2 | surviving_replicate3 | dead_replicate1 | dead_replicate2 | dead_replicate3 |
|---------|----------------------|----------------------|----------------------|-----------------|-----------------|-----------------|
| 126     | Control              | Control              | Control              | Control         | Control         | Control         |
| 127N    | 8-azaguanine         | 8-azaguanine         | 8-azaguanine         | 8-azaguanine    | 8-azaguanine    | 8-azaguanine    |
| 127C    | Raltitrexed          | Raltitrexed          | Raltitrexed          | Raltitrexed     | Raltitrexed     | Raltitrexed     |
| 128N    | Topotecan            | Topotecan            | Topotecan            | Topotecan       | Topotecan       | Topotecan       |
| 128C    | Floxuridine          | Floxuridine          | Floxuridine          | Floxuridine     | Floxuridine     | Floxuridine     |
| 129N    | Nutlin               | Nutlin               | Nutlin               | Nutlin          | Nutlin          | Nutlin          |
| 129C    | Dasatanib            | Dasatanib            | Dasatanib            | Dasatanib       | Dasatanib       | Dasatanib       |
| 130N    | Gefitinib            | Gefitinib            | Gefitinib            | Gefitinib       | Gefitinib       | Gefitinib       |
| 130C    | Vincristine          | Vincristine          | Vincristine          | Vincristine     | Vincristine     | Vincristine     |
| 131     | Bortezomib           | Bortezomib           | Bortezomib           | Bortezomib      | Bortezomib      | Bortezomib      |

## 2021-03-02 Hello Life and Death again...

Starting up a .md log instead of .rmd log.


Links for understanding TMT:

https://www.researchgate.net/figure/Principle-of-multiplexed-proteomics-A-Proteins-from-multiple-conditions-replicates_fig1_327239929

https://idearesourceproteomics.org/wp-content/uploads/2017/03/Tandem-Mass-Tag-TMT-Labeling-Workflow.pdf (p. 5)

So, all the proteins in a sample get a TMT-tag. The MS1 scan contains multiplexed peaks, performing DDA on a peak gives us a MS2 with TMT-tag peaks. Since the TMT-tag have different mass they will tell us the quantity of the peptides for each sample. The MS2 peaks are the same for each sample.





