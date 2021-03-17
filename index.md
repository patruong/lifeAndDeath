# About

## Relevant Papers
[ProTargetMiner - Saei et al. 2019](https://www.nature.com/articles/s41467-019-13582-8)
[Comparative Proteomics - Saei et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29572246/)



# Blog

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





