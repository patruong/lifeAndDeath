if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")

library(ropls)
rm(list=ls())

data(sacurine)
names(sacurine)

attach(sacurine)

#view(dataMatrix)

#heatmap(dataMatrix)
datMatrix <- sacurine$dataMatrix

# PCA
sacurine.pca <- opls(datMatrix)
sacurine.pca <- opls(dataMatrix)

genderFc <- sampleMetadata[, "gender"]
plot(sacurine.pca,
     typeVc = "x-score",
     parAsColFcVn = genderFc)

plot(sacurine.pca,
     typeVc = "x-score",
     parAsColFcVn = genderFc,
     parLabVc = as.character(sampleMetadata[, "age"]),
     parPaletteVc = c("green4", "magenta"))


# PLS-DA
sacurine.plsda <- opls(dataMatrix, genderFc)

# OPLS-DA
sacurine.oplsda <- opls(dataMatrix, genderFc,
                        predI = 1, orthoI = NA,
                        subset = "odd")

trainVi <- getSubsetVi(sacurine.oplsda)
table(genderFc[trainVi], fitted(sacurine.oplsda))

table(genderFc[-trainVi],
      predict(sacurine.oplsda, dataMatrix[-trainVi, ]))





