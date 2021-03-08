# About

## Relevant Papers
[ProTargetMiner - Saei et al. 2019](https://www.nature.com/articles/s41467-019-13582-8)
[Comparative Proteomics - Saei et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29572246/)



# Blog

## 2021-03-03 Table of experiment


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




